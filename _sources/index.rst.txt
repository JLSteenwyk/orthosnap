.. image:: _static/img/logo_no_subtext.png
   :width: 55%
   :align: center
   :target: https://jlsteenwyk.com/orthosnap

^^^^^

OrthoSNAP is a command-line tool for increasing the size of molecular evolution datasets that
can then be used for diverse studies including phylogenomics and genome-wide surveys of positive
selection.

If you found OrthoSNAP useful, please cite *OrthoSNAP: a tree splitting and pruning algorithm
for retrieving single-copy orthologs from gene family trees*. bioRxiv. doi: |doiLink|_.

.. _doiLink: https://www.biorxiv.org/content/10.1101/2021.10.30.466607v1
.. |doiLink| replace:: 10.1101/2021.10.30.466607

|

Quick Start
-----------
**1) Installation**

To install using *pip*, we strongly recommend building a virtual environment to avoid 
software dependency issues. To do so, execute the following commands:

.. code-block:: shell

	# create virtual environment
	python -m venv .venv
	# activate virtual environment
	source .venv/bin/activate
	# install orthosnap
	pip install orthosnap

**Note, the virtual environment must be activated to use orthosnap.**

After using orthosnap, you may wish to deactivate your virtual environment and can do so using the following command:

.. code-block:: shell

	# deactivate virtual environment
	deactivate

|

To install via anaconda, execute the following command:

.. code-block:: shell

	conda install -c jlsteenwyk orthosnap

Visit here for more information:
https://anaconda.org/JLSteenwyk/orthosnap

|

Similarly, to install from source, we strongly recommend using a virtual environment. To do so, use the
following commands:

.. code-block:: shell

	# download
	git clone https://github.com/JLSteenwyk/orthosnap.git
	cd orthosnap/
	# create virtual environment
	python -m venv .venv
	# activate virtual environment
	source .venv/bin/activate
	# install
	make install

To deactivate your virtual environment, use the following command:

.. code-block:: shell

	# deactivate virtual environment
	deactivate

**Note, the virtual environment must be activated to use orthosnap.**

|

**2) Usage**

Get the help message from orthosnap:

.. code-block:: shell

	orthosnap -h

|

^^^^

.. toctree::
	:maxdepth: 4

	about/index
	usage/index
	tutorial/index
	change_log/index
	other_software/index
	faq/index

^^^^

