========
psap
========


.. image:: https://img.shields.io/pypi/v/psap_cli.svg
        :target: https://pypi.python.org/pypi/psap_cli

.. image:: https://img.shields.io/travis/tilschaef/psap_cli.svg
        :target: https://travis-ci.com/tilschaef/psap_cli

.. image:: https://readthedocs.org/projects/psap-cli/badge/?version=latest
        :target: https://psap-cli.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

CLI interface for the PSAP classifier, Mierlo, G. van. Predicting protein condensate formation using machine learning (Manuscript in Preparation).


* Free software: MIT license

========
Getting Started
========

1. *Install psap*
--------
.. code-block:: bash
   
   cd psap && python setup.py install
   
2. *Train classifier*
--------
.. code-block:: python

   psap train -f /path/to/peptide-trainingset.fasta  -o /output/directory  
The trained RandomForest classifiers is exported to json format and stored in the output directory.

3. *Predict llps score for peptide instances*
--------
.. code-block:: python

   psap predict -m /path/to/model.json -f /path/to/peptid-testset.fasta -o /output/directory
   
When no model is provided (-m) psap loads the default classifier stored in /data/model.

Optional
-------

*Annotate petides*
--------
.. code-block:: python

   psap annotate -f /path/to/peptide.fasta -o /output/directory    

Annotates a peptide fasta with biochemical features.



Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
