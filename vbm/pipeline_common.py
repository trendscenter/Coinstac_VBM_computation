"""Shared helpers for VBM use-case workflow modules."""

import contextlib
import os


@contextlib.contextmanager
def stdchannel_redirected(stdchannel, dest_filename):
    """Temporarily redirect stdout or stderr."""
    oldstdchannel = None
    dest_file = None
    try:
        oldstdchannel = os.dup(stdchannel.fileno())
        dest_file = open(dest_filename, 'w')
        os.dup2(dest_file.fileno(), stdchannel.fileno())
        yield
    finally:
        if oldstdchannel is not None:
            os.dup2(oldstdchannel, stdchannel.fileno())
        if dest_file is not None:
            dest_file.close()


def create_workflow_input(source, target, source_output, target_input):
    """Collect pipeline nodes and their connection tuple."""
    return (source, target, [(source_output, target_input)])
