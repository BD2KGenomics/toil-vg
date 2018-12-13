#!/usr/bin/env python2.7
"""
Shared stuff between different modules in this package.  Some
may eventually move to or be replaced by stuff in toil-lib.
"""
from __future__ import print_function
import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import json, timeit, errno
import fcntl
import select
import time
import threading
from uuid import uuid4
import pkg_resources, tempfile, datetime
import logging
from distutils.spawn import find_executable
import collections

from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger
from toil.lib.docker import dockerCall, dockerCheckOutput, apiDockerCall
from toil_vg.singularity import singularityCall, singularityCheckOutput
from toil_vg.iostore import IOStore

logger = logging.getLogger(__name__)

# We need to fiddle with os.environ and we don't want to step on ourselves
environment_lock = threading.Lock()

def test_docker():
    """
    Return true if Docker is available on this machine, and False otherwise.
    """
    
    # We don't actually want any Docker output.
    nowhere = open(os.devnull, 'wb')
    
    try:
        # Run Docker
        # TODO: implement around dockerCall somehow?
        subprocess.check_call(['docker', 'version'], stdout=nowhere, stderr=nowhere)
        # And report that it worked
        return True
    except:
        # It didn't work, so we can't use Docker
        return False

def add_container_tool_parse_args(parser):
    """ centralize shared container options and their defaults """

    parser.add_argument("--vg_docker", type=str,
                        help="Docker image to use for vg")
    parser.add_argument("--container", default=None, choices=['Docker', 'Singularity', 'None'],
                       help="Container type used for running commands. Use None to "
                       " run locally on command line")    

def add_common_vg_parse_args(parser):
    """ centralize some shared io functions and their defaults """
    parser.add_argument('--config', default=None, type=str,
                        help='Config file.  Use toil-vg generate-config to see defaults/create new file')
    parser.add_argument('--whole_genome_config', action='store_true',
                        help='Use the default whole-genome config (as generated by toil-vg config --whole_genome)')
    parser.add_argument("--force_outstore", action="store_true",
                        help="use output store instead of toil for all intermediate files (use only for debugging)")
    parser.add_argument("--realTimeStderr", action="store_true",
                        help="print stderr from all commands through the realtime logger")                        
    
def get_container_tool_map(options):
    """ convenience function to parse the above _container options into a dictionary """

    cmap = [dict(), options.container]
    cmap[0]["vg"] = options.vg_docker
    cmap[0]["bcftools"] = options.bcftools_docker
    cmap[0]["tabix"] = options.tabix_docker
    cmap[0]["bgzip"] = options.tabix_docker
    cmap[0]["jq"] = options.jq_docker
    cmap[0]["rtg"] = options.rtg_docker
    cmap[0]["pigz"] = options.pigz_docker
    cmap[0]["samtools"] = options.samtools_docker
    cmap[0]["bwa"] = options.bwa_docker
    cmap[0]["Rscript"] = options.r_docker
    cmap[0]["vcfremovesamples"] = options.vcflib_docker
    cmap[0]["freebayes"] = options.freebayes_docker
    cmap[0]["Platypus.py"] = options.platypus_docker
    cmap[0]['hap.py'] = options.happy_docker
    cmap[0]['bedtools'] = options.bedtools_docker
    cmap[0]['R'] = options.sveval_docker
    cmap[0]['bedops'] = options.bedops_docker
     
    # to do: could be a good place to do an existence check on these tools

    return cmap

def toil_call(job, context, cmd, work_dir, out_path = None, out_append = False):
    """ use to run a one-job toil workflow just to call a command
    using context.runner """
    if out_path:
        open_flag = 'a' if out_append is True else 'w'
        with open(os.path.abspath(out_path), open_flag) as out_file:
            context.runner.call(job, cmd, work_dir=work_dir, outfile=out_file)
    else:
        context.runner.call(job, cmd, work_dir=work_dir)        
    
class ContainerRunner(object):
    """ Helper class to centralize container calling.  So we can toggle both
Docker and Singularity on and off in just one place.
to do: Should go somewhere more central """
    def __init__(self, container_tool_map = [{}, None], realtime_stderr=False):
        # this maps a command to its full docker name
        #   the first index is a dictionary containing docker tool names
        #   the second index is a string that represents which container
        #   support to use.
        # example:  docker_tool_map['vg'] = 'quay.io/ucsc_cgl/vg:latest'
        #           container_support = 'Docker'
        self.docker_tool_map = container_tool_map[0]
        self.container_support = container_tool_map[1]
        self.realtime_stderr = realtime_stderr

    def container_for_tool(self, name):
        """
        Return Docker, Singularity or None, which is how call() would be run
        on the given tool
        """
        if self.container_support == 'Docker' and name in self.docker_tool_map and\
           self.docker_tool_map[name] and self.docker_tool_map[name].lower() != 'none':
            return 'Docker'
        elif self.container_support == 'Singularity' and name in self.docker_tool_map and\
           self.docker_tool_map[name] and self.docker_tool_map[name].lower() != 'none':
            return 'Singularity'
        else:
            return 'None'

    def call(self, job, args, work_dir = '.' , outfile = None, errfile = None,
             check_output = False, tool_name=None):
        """
        
        Run a command. Decide to use a container based on whether the tool
        (either the tool of the first command, or the tool named by tool_name)
        its in the container engine's tool map. Can handle args being either a
        single list of string arguments (starting with the binary name) or a
        list of such lists (in which case a pipeline is run, using no more than
        one container).
        
        Redirects standard output and standard error to the file objects outfile
        and errfile, if specified. If check_output is true, the call will block,
        raise an exception on a nonzero exit status, and return standard
        output's contents.
        
        """
        # from here on, we assume our args is a list of lists
        if len(args) == 0 or len(args) > 0 and type(args[0]) is not list:
            args = [args]
        # convert everything to string
        for i in range(len(args)):
            args[i] = [str(x) for x in args[i]]
        name = tool_name if tool_name is not None else args[0][0]

        # optionally log stderr to the realtime logger by making a pipe and
        # logging the output in a forked process
        # todo: duplicate errfile to logger rather than ignoring when errfile not None
        if self.realtime_stderr and not errfile:
            # Make our pipe
            rfd, wfd = os.pipe()
            rfile = os.fdopen(rfd, 'r', 0)
            wfile = os.fdopen(wfd, 'w', 0)
            # Fork our child process (pid == 0) to catch stderr and log it
            pid = os.fork()
            if pid == 0:
                wfile.close()
                while 1:
                    data = rfile.readline()
                    if not data:
                        break
                    RealtimeLogger.info('(stderr) {}'.format(data.strip()))
                os._exit(0)
            else:
                assert pid > 0
                # main process carries on, but sending stderr to the pipe
                rfile.close()
                # note that only call_directly below actually does anything with errfile at the moment
                errfile = wfile

        container_type = self.container_for_tool(name)
                
        if container_type == 'Docker':
            return self.call_with_docker(job, args, work_dir, outfile, errfile, check_output, tool_name)
        elif container_type == 'Singularity':
            return self.call_with_singularity(job, args, work_dir, outfile, errfile, check_output, tool_name)
        else:
            return self.call_directly(args, work_dir, outfile, errfile, check_output)
        
    def call_with_docker(self, job, args, work_dir, outfile, errfile, check_output, tool_name): 
        """
        
        Thin wrapper for docker_call that will use internal lookup to
        figure out the location of the docker file.  Only exposes docker_call
        parameters used so far.  expect args as list of lists.  if (toplevel)
        list has size > 1, then piping interface used
        
        Does support redirecting output to outfile, unless check_output is
        used, in which case output is captured.
        
        """

        RealtimeLogger.info("Docker Run: {}".format(" | ".join(" ".join(x) for x in args)))
        start_time = timeit.default_timer()

        # we use the first argument to look up the tool in the docker map
        # but allow overriding of this with the tool_name parameter
        name = tool_name if tool_name is not None else args[0][0]
        tool = self.docker_tool_map[name]

        # We keep an environment dict
        environment = {}
        
        # And an entry point override
        entrypoint = None
        
        # And a volumes dict for mounting
        volumes = {}
        
        # And a working directory override
        working_dir = None

        # breaks Rscript.  Todo: investigate how general this actually is
        if name != 'Rscript':
            # vg uses TMPDIR for temporary files
            # this is particularly important for gcsa, which makes massive files.
            # we will default to keeping these in our working directory
            environment['TMPDIR'] = '.'
            
        if name == 'Rscript':
            # The R dockers by default want to install packages in non-writable directories. Sometimes.
            # Make sure a writable directory which exists is used.
            environment['R_LIBS']='/tmp'

        if name == 'vg':
            environment['VG_FULL_TRACEBACK'] = '1'

        # ugly hack for platypus, as default container doesn't have executable in path
        if tool == 'quay.io/biocontainers/platypus-variant:0.8.1.1--htslib1.7_1' and \
           args[0][0] == 'Platypus.py':
            args[0][0] = '/usr/local/share/platypus-variant-0.8.1.1-1/Platypus.py'

        # Force all dockers to run sort in a consistent way
        environment['LC_ALL'] = 'C'

        # set our working directory map
        if work_dir is not None:
            volumes[os.path.abspath(work_dir)] = {'bind': '/data', 'mode': 'rw'}
            working_dir = '/data'
            
        if outfile is not None:
            # We need to send output to a file object

            assert(not check_output)
            
            # We can't just redirect stdout of the container from the API, so
            # we do something more complicated.
            
            # Set up a FIFO to receive it
            fifo_dir = tempfile.mkdtemp()
            fifo_host_path = os.path.join(fifo_dir, 'stdout.fifo')
            os.mkfifo(fifo_host_path)
            
            # Mount the FIFO in the container.
            # The container doesn't actually have to have the mountpoint directory in its filesystem.
            volumes[fifo_dir] = {'bind': '/control', 'mode': 'rw'}
            
            # Redirect the command output by tacking on another pipeline stage
            parameters = args + [['dd', 'of=/control/stdout.fifo']]
            
            
            # Start the container detached so we don't wait on it
            container = apiDockerCall(job, tool, parameters,
                                      volumes=volumes,
                                      working_dir=working_dir,
                                      entrypoint=entrypoint,
                                      environment=environment,
                                      detach=True)
            
            RealtimeLogger.debug("Asked for container {}".format(container.id))
            
            # If the Docker container goes badly enough, it may not even
            # open the other end of the FIFO. So we can't just wait for it
            # to EOF before checking on the Docker.
            
            # Open the FIFO into nonblocking mode. See
            # <https://stackoverflow.com/a/5749687> and
            # <http://shallowsky.com/blog/programming/python-read-characters.html>
            fifo_fd = os.open(fifo_host_path, os.O_RDONLY | os.O_NONBLOCK)
            
            # Now read ought to throw if there is no data. But
            # <https://stackoverflow.com/q/38843278> and some testing suggest
            # that this doesn't happen, and it just looks like EOF. So we will
            # watch out for that.
            
            try:
                # Prevent leaking FDs
            
                # If this is set, and there is no data in the pipe, decide that no data is coming
                last_chance = False
                # If this is set, we have seen data in the pipe, so the other
                # end must have opened it and will eventually close it if it
                # doesn't run forever.
                saw_data = False
                
                while True:
                    # While there still might be data in the pipe
                    
                    # Select on the pipe with a timeout, so we don't spin constantly waiting for data
                    can_read, can_write, had_error = select.select([fifo_fd], [], [fifo_fd], 10)
                    
                    if len(can_read) > 0 or len(had_error) > 0:
                        # There is data available or something else weird about our FIFO.
                    
                        try:
                            # Do a nonblocking read. Since we checked with select we never should get "" unless there's an EOF.
                            data = os.read(fifo_fd, 4096)
                            
                            if data == "":
                                # We didn't throw and we got nothing, so it must be EOF.
                                RealtimeLogger.debug("Got EOF")
                                break
                                
                        except OSError as err:
                            if err.errno in [errno.EAGAIN, errno.EWOULDBLOCK]:
                                # There is no data right now
                                data = None
                            else:
                                # Something else has gone wrong
                                raise err
                                
                    else:
                        # There is no data available. Don't even try to read. Treat it as if a read refused to block.
                        data = None
                    
                    
                        
                    if data is not None:
                        # Send our data to the outfile
                        outfile.write(data)
                        saw_data = True
                    elif not saw_data:
                        # We timed out and there has never been any data. Maybe the container has died/never started?
                        
                        if last_chance:
                            # The container has been dead for a while and nothing has arrived yet. Assume no data is coming.
                            RealtimeLogger.warning("Giving up on output form container {}".format(container.id))
                            break
                        
                        # Otherwise, check on it
                        container.reload()
                        
                        if container.status not in ['created', 'restarting', 'running', 'removing']:
                            # The container has stopped. So what are we doing waiting around for it?
                            
                            # Wait one last time for any lingering data to percolate through the FIFO
                            time.sleep(10)
                            last_chance = True
                            continue
                            
            finally:
                # No matter what happens, close our end of the FIFO
                os.close(fifo_fd)
                    
                        
            # Now our data is all sent.
            # Wait on the container and get its return code.
            return_code = container.wait()
            
            os.unlink(fifo_host_path)
            os.rmdir(fifo_dir)
            
        else:
            # No piping needed.
        
            if len(args) == 1:
                # split off first argument as entrypoint (so we can be oblivious as to whether
                # that happens by default)
                parameters = [] if len(args[0]) == 1 else args[0][1:]
                entrypoint = args[0][0]
            else:
                # can leave as is for piped interface which takes list of args lists
                # and doesn't worry about entrypoints since everything goes through bash -c
                # todo: check we have a bash entrypoint!
                parameters = args
                
                
            # Run the container and dump the logs if it fails.
            container = apiDockerCall(job, tool, parameters,
                                      volumes=volumes,
                                      working_dir=working_dir,
                                      entrypoint=entrypoint,
                                      environment=environment,
                                      detach=True)
                                      
            # Wait on the container and get its return code.
            return_code = container.wait()
            
        # When we get here, the container has been run, and stdout is either in the file object we sent it to or in the Docker logs.
        # stderr is always in the Docker logs.
            
        if return_code != 0:
            # What were we doing?
            command = " | ".join(" ".join(x) for x in args)
        
            # Dump logs
            RealtimeLogger.error("Docker container for command {} failed with code {}".format(command, return_code))
            RealtimeLogger.error("Dumping stderr...")
            for line in container.logs(stderr=True, stdout=False, stream=True):
                # Trim trailing \n
                RealtimeLogger.error(line[:-1])
                
            if not check_output and outfile is None:
                # Dump stdout as well, since it's not something the caller wanted as data
                RealtimeLogger.error("Dumping stdout...")
                for line in container.logs(stderr=False, stdout=True, stream=True):
                    # Trim trailing \n
                    RealtimeLogger.error(line[:-1])
        
            # Raise an error if it's not sucess
            raise RuntimeError("Docker container for command {} failed with code {}".format(command, return_code))
        elif errfile:
            # user wants stderr even if no crash
            for line in container.logs(stderr=True, stdout=False, stream=True):
                errfile.write(line)
        
        if check_output:
            # We need to collect the output. We grab it from Docker's handy on-disk buffer.
            # TODO: Bad Things can happen if the container logs too much.
            captured_stdout = container.logs(stderr=False, stdout=True)
            
        end_time = timeit.default_timer()
        run_time = end_time - start_time
        RealtimeLogger.info("Successfully docker ran {} in {} seconds.".format(
            " | ".join(" ".join(x) for x in args), run_time))
        
        if outfile:
            outfile.flush()
            os.fsync(outfile.fileno())

        if check_output is True:
            return captured_stdout
    
    def call_with_singularity(self, job, args, work_dir, outfile, errfile, check_output, tool_name): 
        """ Thin wrapper for singularity_call that will use internal lookup to
        figure out the location of the singularity file.  Only exposes singularity_call
        parameters used so far.  expect args as list of lists.  if (toplevel)
        list has size > 1, then piping interface used """

        RealtimeLogger.info("Singularity Run: {}".format(" | ".join(" ".join(x) for x in args)))
        start_time = timeit.default_timer()

        # we use the first argument to look up the tool in the singularity map
        # but allow overriding of this with the tool_name parameter
        name = tool_name if tool_name is not None else args[0][0]
        tool = self.docker_tool_map[name]

        parameters = args[0] if len(args) == 1 else args
        
        # Get a lock on the environment
        global environment_lock
        with environment_lock:
            # TODO: We can't stop other threads using os.environ or subprocess or w/e on their own

            # Set the locale to C for consistent sorting, and activate vg traceback     
            update_env = {'LC_ALL' : 'C', 'VG_FULL_TRACEBACK': '1', 'TMPDIR': '.'}
            old_env = {}
            for env_name, env_val in update_env.items():
                old_env[env_name] = os.environ.get(env_name)
                os.environ[env_name] = env_val
            
            if check_output is True:
                ret = singularityCheckOutput(job, tool, parameters=parameters, workDir=work_dir)
            else:
                ret = singularityCall(job, tool, parameters=parameters, workDir=work_dir, outfile = outfile)
            
            # Restore old locale and vg traceback
            for env_name, env_val in update_env.items():
                if old_env[env_name] is not None:
                    os.environ[env_name] = old_env[env_name]
                else:
                    del os.environ[env_name]
        
        end_time = timeit.default_timer()
        run_time = end_time - start_time
        RealtimeLogger.info("Successfully singularity ran {} in {} seconds.".format(
            " | ".join(" ".join(x) for x in args), run_time))

        if outfile:
            outfile.flush()
            os.fsync(outfile.fileno())
        
        return ret

    def call_directly(self, args, work_dir, outfile, errfile, check_output):
        """ Just run the command without docker """

        RealtimeLogger.info("Run: {}".format(" | ".join(" ".join(x) for x in args)))
        start_time = timeit.default_timer()

        # Set up the child's environment
        global environment_lock
        with environment_lock:
            my_env = os.environ.copy()
            
        # vg uses TMPDIR for temporary files
        # this is particularly important for gcsa, which makes massive files.
        # we will default to keeping these in our working directory
        my_env['TMPDIR'] = '.'
        
        # Set the locale to C for consistent sorting
        my_env['LC_ALL'] = 'C'

        # Set the vg traceback
        my_env['VG_FULL_TRACEBACK'] = '1'

        procs = []
        for i in range(len(args)):
            stdin = procs[i-1].stdout if i > 0 else None
            if i == len(args) - 1 and outfile is not None:
                stdout = outfile
            else:
                stdout = subprocess.PIPE

            try:
                procs.append(subprocess.Popen(args[i], stdout=stdout, stderr=errfile,
                                              stdin=stdin, cwd=work_dir, env=my_env))
            except OSError as e:
                # the default message: OSError: [Errno 13] Permission denied is a bit cryptic
                # so we print something a bit more explicit if a command isn't found
                if e.errno in [2,13] and not find_executable(args[i][0]):
                    raise RuntimeError('Command not found: {}'.format(args[i][0]))
                else:
                    raise e
            
        for p in procs[:-1]:
            p.stdout.close()

        output, errors = procs[-1].communicate()
        for i, proc in enumerate(procs):
            sts = proc.wait()
            if sts != 0:
                raise Exception("Command {} returned with non-zero exit status {}".format(
                    " ".join(args[i]), sts))

        end_time = timeit.default_timer()
        run_time = end_time - start_time
        RealtimeLogger.info("Successfully ran {} in {} seconds.".format(
            " | ".join(" ".join(x) for x in args), run_time))            
        
        if outfile:
            outfile.flush()
            os.fsync(outfile.fileno())

        if check_output:
            return output

def get_vg_script(job, runner, script_name, work_dir):
    """
    getting the path to a script in vg/scripts is different depending on if we're
    in docker or not.  wrap logic up here, where we get the script from wherever it
    is then put it in the work_dir
    """
    vg_container_type = runner.container_for_tool('vg')

    if vg_container_type != 'None':
        # we copy the scripts out of the container, assuming vg is at /vg
        cmd = ['cp', os.path.join('/vg', 'scripts', script_name), '.']
        runner.call(job, cmd, work_dir = work_dir, tool_name='vg')
    else:
        # we copy the script from the vg directory in our PATH
        scripts_path = os.path.join(os.path.dirname(find_executable('vg')), '..', 'scripts')
        shutil.copy2(os.path.join(scripts_path, script_name), os.path.join(work_dir, script_name))
    return os.path.join(work_dir, script_name)
                            
def get_files_by_file_size(dirname, reverse=False):
    """ Return list of file paths in directory sorted by file size """

    # Get list of files
    filepaths = []
    for basename in os.listdir(dirname):
        filename = os.path.join(dirname, basename)
        if os.path.isfile(filename):
            filepaths.append(filename)

    # Re-populate list with filename, size tuples
    for i in xrange(len(filepaths)):
        filepaths[i] = (filepaths[i], os.path.getsize(filepaths[i]))

    return filepaths

def make_url(path):
    """ Turn filenames into URLs, whileleaving existing URLs alone """
    # local path
    if ':' not in path:
        return 'file://' + os.path.abspath(path)
    else:
        return path
    
def require(expression, message):
    if not expression:
        raise Exception('\n\n' + message + '\n\n')

def parse_id_ranges(job, id_ranges_file_id):
    """Returns list of triples chrom, start, end
    """
    work_dir = job.fileStore.getLocalTempDir()
    id_range_file = os.path.join(work_dir, 'id_ranges.tsv')
    job.fileStore.readGlobalFile(id_ranges_file_id, id_range_file)
    return parse_id_ranges_file(id_range_file)

def parse_id_ranges_file(id_ranges_filename):
    """Returns list of triples chrom, start, end
    """
    id_ranges = []
    with open(id_ranges_filename) as f:
        for line in f:
            toks = line.split()
            if len(toks) == 3:
                id_ranges.append((toks[0], int(toks[1]), int(toks[2])))
    return id_ranges

def remove_ext(string, ext):
    """
    Strip a suffix from a string. Case insensitive.
    """
    # See <https://stackoverflow.com/a/18723694>
    if string.lower().endswith(ext.lower()):
        return string[:-len(ext)]
    else:
        return string

class TimeTracker:
    """ helper dictionary to keep tabs on several named runtimes. """
    def __init__(self, name = None):
        """ create. optionally start a timer"""
        self.times = collections.defaultdict(float)
        self.running = {}
        if name:
            self.start(name)
    def start(self, name):
        """ start a timer """
        assert name not in self.running
        self.running[name] = timeit.default_timer()
    def stop(self, name = None):
        """ stop a timer. if no name, do all running """
        names = [name] if name else self.running.keys()
        ti = timeit.default_timer()
        for name in names:
            self.times[name] += ti - self.running[name]
            del self.running[name]
    def add(self, time_dict):
        """ add in all times from another TimeTracker """
        for key, value in time_dict.times.items():
            self.times[key] += value
    def total(self, names = None):
        if not names:
            names = self.times.keys()
        return sum([self.times[name] for name in names])
    def names(self):
        return self.times.keys()
        
        
def run_concat_lists(job, *args):
    """
    Toil job to join all the given lists and return the merged list.
    """
    
    concat = []
    for input_list in args:
        concat += input_list
    return concat
    
def parse_plot_set(plot_set_string):
    """
    
    Given one of the string arguments to the --plot-sets option, parse out a
    data structure representing which conditions ought to be compared against
    each other, and what those comparison plots/tables should be called.

    The syntax of a plot set is [title:]condition[,condition[,condition...]].
    
    The first condition is the comparison baseline, when applicable.
    
    Returns a tuple of a plot set title, or None if unspecified, and a list of
    condition names.
    
    """
    
    colon_pos = plot_set_string.find(':')
    
    if colon_pos != -1:
        # Pull out the title before the colon
        title = plot_set_string[0:colon_pos]
        # And the rest of the specifier after it
        plot_set_string = plot_set_string[colon_pos + 1:]
    else:
        # No title given
        title = None
        
    # Return the title and condition list tuple
    return (title, plot_set_string.split(','))
        
    
def parse_plot_sets(plot_sets_list):
    """
    
    Given a list of plot set strings, parses each with parse_plot_set.
    
    Returns a list of tuples. Each tuple is a plot set title, or None if no
    title is to be applied, and a list of condition names, or None if all
    conditions are to be included.
    
    If no plot sets are specified in the list, returns a single plot set for
    all conditions.
    
    """
    
    plot_sets = [parse_plot_set(spec) for spec in plot_sets_list]
    if len(plot_sets) == 0:
        # We want to plot everything together
        # We use the special None value instead of a condition list to request that.
        # And we use None for the title.
        plot_sets = [(None, None)]
    
    return plot_sets
    
def title_to_filename(kind, i, title, extension):
    """
    Given the kind of thing you want to save ('table', 'plot-qq', etc.), the
    number of the thign out of all things of that type, the human-readable
    title ('All Conditions vs. Whatever'), and an extansion, come up with a
    safe filename to save the plot under.
    
    The title may be None, in which case it is ommitted.
    
    The extension may be None, in which case it is omitted.
    """
    
    if title is not None:
        # Filter down to good filename characters as in https://stackoverflow.com/a/7406369
        safe_title = ''.join((c for c in title if c.isalnum()))
    else:
        safe_title = None
        
    # The name always includes the kind of thing
    part_list = [kind]
    
    if i != 0:
        # Include number padded to 2 digits only if nonzero. If not included,
        # the next dash will sort before the other numbers. TODO: Add the 00 in
        # here and break backward compatibility with everything looking for the
        # output files.
        part_list.append('-{:02d}'.format(i))
        
    if title is not None:
        # Filter down to good filename characters as in https://stackoverflow.com/a/7406369
        safe_title = ''.join((c for c in title if c.isalnum()))
        part_list.append('-{}'.format(safe_title))
        
    if extension is not None:
        part_list.append('.{}'.format(extension))
        
    return ''.join(part_list)

