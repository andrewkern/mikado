Installation
============

Requirements
------------

- Python 3.12 or later
- NumPy
- SciPy
- Typer
- Rich

Installing with pip
-------------------

The easiest way to install MKado is via pip:

.. code-block:: bash

   pip install mkado

Installing with uv (Recommended)
--------------------------------

For faster installation and dependency management, use `uv <https://docs.astral.sh/uv/>`_:

.. code-block:: bash

   uv pip install mkado

Development Installation
------------------------

To install for development:

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/andrewkern/mkado.git
   cd mkado

   # Install with uv (recommended)
   uv sync

   # Or install with pip in editable mode
   pip install -e .

Running Tests
-------------

After development installation:

.. code-block:: bash

   # Run all tests
   uv run pytest

   # Run with coverage
   uv run pytest --cov=mkado

   # Run specific test file
   uv run pytest tests/test_mk_test.py -v

Verifying Installation
----------------------

After installation, verify it works:

.. code-block:: bash

   mkado --help

You should see the help output listing available commands.
