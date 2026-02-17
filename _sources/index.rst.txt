.. image:: _static/img/logo_no_subtext.png
   :width: 55%
   :align: center
   :target: https://jlsteenwyk.com/orthosnap

OrthoSNAP is a command-line tool for retrieving single-copy orthologous subgroups (SNAP-OGs)
from larger gene families.

If you found OrthoSNAP useful, please cite *OrthoSNAP: a tree splitting and pruning algorithm
for retrieving single-copy orthologs from gene family trees*. PLOS Biology. doi: |doiLink|_.

.. _doiLink: https://jlsteenwyk.com/publication_pdfs/2022_Steenwyk_etal_PLoS_Biology.pdf
.. |doiLink| replace:: 10.1371/journal.pbio.3001827

Quick Start
-----------

1. Install with pip
^^^^^^^^^^^^^^^^^^^

We recommend using a virtual environment.

.. code-block:: shell

   python -m venv .venv
   source .venv/bin/activate
   pip install orthosnap

Deactivate when finished:

.. code-block:: shell

   deactivate

2. Install with conda
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

   conda install -c jlsteenwyk orthosnap

More details: https://anaconda.org/JLSteenwyk/orthosnap

3. Install from source
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

   git clone https://github.com/JLSteenwyk/orthosnap.git
   cd orthosnap
   python -m venv .venv
   source .venv/bin/activate
   make install

4. Show CLI help
^^^^^^^^^^^^^^^^

.. code-block:: shell

   orthosnap -h

.. toctree::
   :maxdepth: 4

   about/index
   usage/index
   tutorial/index
   change_log/index
   other_software/index
   faq/index
