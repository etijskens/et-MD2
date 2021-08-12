# -*- coding: utf-8 -*-

"""
Package et_md2
==============

Top-level package for et_md2.
"""

__version__ = "0.0.0"

try:
    import et_md2.potentials.cpp
except ModuleNotFoundError as e:
    # Try to build this binary extension:
    from pathlib import Path
    import click
    from et_micc2.project import auto_build_binary_extension
    msg = auto_build_binary_extension(Path(__file__).parent, 'potentials/cpp')
    if not msg:
        import et_md2.potentials.cpp
    else:
        click.secho(msg, fg='bright_red')

import et_md2.potentials

try:
    import et_md2.interactions.cpp
except ModuleNotFoundError as e:
    # Try to build this binary extension:
    from pathlib import Path
    import click
    from et_micc2.project import auto_build_binary_extension
    msg = auto_build_binary_extension(Path(__file__).parent, 'interactions/cpp')
    if not msg:
        import et_md2.interactions.cpp
    else:
        click.secho(msg, fg='bright_red')

import et_md2.interactions.lj

import et_md2.interactions

try:
    import et_md2.verletlist.c_vl
except ModuleNotFoundError as e:
    # Try to build this binary extension:
    from pathlib import Path
    import click
    from et_micc2.project import auto_build_binary_extension
    msg = auto_build_binary_extension(Path(__file__).parent, 'verletlist/impl_cpp')
    if not msg:
        import et_md2.verletlist.c_vl
    else:
        click.secho(msg, fg='bright_red')


import et_md2.verletlist

import et_md2.atoms

