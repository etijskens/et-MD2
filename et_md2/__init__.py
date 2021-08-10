# -*- coding: utf-8 -*-

"""
Package et_md2
==============

Top-level package for et_md2.
"""

__version__ = "0.0.0"

try:
    import et_md2.verletlist.impl_cpp
except ModuleNotFoundError as e:
    # Try to build this binary extension:
    from pathlib import Path
    import click
    from et_micc2.project import auto_build_binary_extension
    msg = auto_build_binary_extension(Path(__file__).parent, 'verletlist/impl_cpp')
    if not msg:
        import et_md2.verletlist.impl_cpp
    else:
        click.secho(msg, fg='bright_red')


import et_md2.verletlist

import et_md2.atoms


# Your code here...