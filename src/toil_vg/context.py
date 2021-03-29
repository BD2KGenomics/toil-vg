#!/usr/bin/env python2.7
"""
context.py: Defines a toil_vg context, which contains (and hides) all the config
file values, IOStore dumping parameters, and other things that need to be passed
around but that toil-vg users shouldn't generally need to dinker with.

Instead of dropping the container runner into the command-line options
namespace, we keep them both in here.

"""

import tempfile
import datetime
import pkg_resources
import copy
import os
import os.path

from argparse import Namespace

from toil.common import Toil
from toil.job import Job

from toil_vg.vg_config import apply_config_file_args
from toil_vg.vg_common import ContainerRunner, get_container_tool_map
from toil_vg.iostore import IOStore

def run_write_info_to_outstore(job, context, argv):
    """ Writing to the output is still problematic, especially from within jobs.  So 
    we write some options info into the output store first-thing to trigger errors
    before doing all the compute if possible.  To do this, this job needs to be passed
    to the root of the workflow """
    
    f = tempfile.NamedTemporaryFile(mode='w',delete=True)
    now = datetime.datetime.now()
    if argv:
        f.write('{}\n\n'.format(' '.join(argv)))
    try:
        version = pkg_resources.get_distribution('toil-vg').version
    except pkg_resources.DistributionNotFound:
        # TODO: pkg_resources gets upset by Kubernetes Toil hot deployment (but not Mesos hot deployment?)
        version = 'unknown'
    f.write('{}\ntoil-vg version {}\nConfiguration:\n'.format(now,
        version))

    for key, val in list(context.config.__dict__.items()):
        f.write('{}: {}\n'.format(key, val))
    f.flush()
    context.get_out_store().write_output_file(f.name, 'toil-vg-{}.txt'.format(
        argv[1] if argv and len(argv) > 1 else 'info'))
    f.close()

class Context(object):
    """
    Represents a toil-vg context, necessary to use the library.
    """
    
    def __init__(self, out_store=None, overrides=Namespace()):
        """
        Make a new context, so we can run the library.
        
        Takes an optional Namespace of overrides for default toil-vg and Toil
        configuration values.
        
        Overrides can also have a "config" key, in which case that config file
        will be loaded.
        
        """
        
        # Load configuration and apply overrides. If the overrides are from
        # command-line options, we might also get a bunch of tool-specific
        # fields.
        self.config = apply_config_file_args(overrides)
        
        # Make a container runner for running tools
        self.runner = ContainerRunner(container_tool_map=get_container_tool_map(
            self.config), realtime_stderr=self.config.realTimeStderr)
        
        if out_store is not None:
            # We want to dump files to an IOStore
            # Make it an absolute path while we're getting set up.
            self.out_store_string = IOStore.absolute(out_store)            
        else:
            # We don't want to use an out store
            self.out_store_string = None
    
    def get_toil(self, job_store):
        """
        Produce a new Toil object for running Toil jobs, using the configuration
        used when constructing this context, and the given job_store specifier.
        
        Needs to be used in a "with" statement to actually work.
        """
        
        # Get the default Toil options
        toil_options = Job.Runner.getDefaultOptions(job_store)
        
        for k, v in list(self.config.__dict__.items()):
            # Blit over all the overrides, some of which will be relevant
            toil_options.__dict__[k] = v
        
        # Make the Toil object and return it.
        # It still needs to be entered as a context manager in order to be used.
        return Toil(toil_options)
            
    def get_out_store(self):
        """
        Return the IOStore to write output files to, or None if they should not
        be written anywhere besides the Toil file store.
        """
        
        if self.out_store_string is not None:
            return IOStore.get(self.out_store_string)
        else:
            return None
            
    def write_intermediate_file(self, job, path, out_store_path = None):
        """
        Write the file at the given path to the given job's Toil FileStore, and
        to the out_store if one is in use and we are trying to dump intermediate
        files.
        
        In the out_store, the file is saved under its basename, in the root
        directory. If out_store_path is set, the file is saved there instead.
        
        Returns the Toil file ID for the written file.
        """
        
        out_store = self.get_out_store()
        if out_store is not None and self.config.force_outstore:
            name = out_store_path if out_store_path else os.path.basename(path)
            # Save to the out_store if it exists
            out_store.write_output_file(path, name)
        
        # Save to Toil
        return job.fileStore.writeGlobalFile(path)
            
    def write_output_file(self, job, path, out_store_path = None):
        """
        
        Write the file at the given path to the given job's Toil FileStore, and
        to the out_store if one is in use.
        
        In the out_store, the file is saved under its basename, in the root
        directory. If out_store_path is set, the file is saved there instead.
        
        Returns the Toil file ID for the written file.
        """
        
        out_store = self.get_out_store()
        if out_store is not None:
            name = out_store_path if out_store_path else os.path.basename(path)
            # Save to the out_store if it exists
            out_store.write_output_file(path, name)
        
        # Save to Toil
        return job.fileStore.writeGlobalFile(path)
            
    def to_options(self, options):
        """
        Turn back into an options-format namespace like toil-vg used to use
        everywhere. Merges with the given real command-line options.
        """
        
        # Copy the options
        options = copy.deepcopy(options)
        
        for key, val in list(self.config.__dict__.items()):
            # Copy over everything from the config to the options namespace.
            # Anything in the config that needed to be overridden by command
            # line options already was overridden in the constructor.
            options.__dict__[key] = val

        # Save the other aspects of the context into the options under their old
        # names.
        options.out_store = self.out_store_string
        options.drunner = self.runner
        options.tool = 'library'
        
        return options
        
        
        
        
        
