"""Example Google style docstrings.

This module demonstrates documentation as specified by the `Google
Python Style Guide`_. Docstrings may extend over multiple lines.
Sections are created with a section header and a colon followed by a
block of indented text.

Example:
    Examples can be given using either the ``Example`` or ``Examples``
    sections. Sections support any reStructuredText formatting, including
    literal blocks::

        $ python example_google.py

Section breaks are created by resuming unindented text. Section breaks
are also implicitly created anytime a new section starts.

Attributes:
    module_level_variable1 (int): Module level variables may be documented in
        either the ``Attributes`` section of the module docstring, or in an
        inline docstring immediately following the variable.

        Either form is acceptable, but the two should not be mixed. Choose
        one convention to document module level variables and be consistent
        with it.

Todo:
    * For module TODOs
    * You have to also use ``sphinx.ext.todo`` extension

.. _Google Python Style Guide:
http://google.github.io/styleguide/pyguide.html

"""
import sys
from pathlib import Path
from flask import Flask, render_template, url_for, request

from FUSION.logger import get_logger
Logger = get_logger("Fusion")


def main():
    from app.server import app
    from FUSION.env import initializer
    from FUSION.arguments import get_parser

    parser = get_parser()
    args = parser.parse_args()

    if Path(sys.argv[0]).name != 'Fusion':
        Logger.setLevel(10)
        for handler in Logger.handlers:
            handler.setLevel(10)

    minknow_host = '127.0.0.1'
    minknow_port = '8000'
    guppy_address = 'ipc:///tmp/.guppy/5555'
    guppy_config = 'dna_r9.4.1_450bps_fast'
    mm_idx = '/home/zhangjy/data/hg38/GRCh38.p13.genome.splice.mmi'
    initializer(minknow_host, minknow_port, guppy_address, guppy_config, mm_idx)

    Logger.info("Everything works fine! Starting dashboard app")
    app.run(debug=False, host="0.0.0.0", port="55280")


if __name__ == '__main__':
    main()
