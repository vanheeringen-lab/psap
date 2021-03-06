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
* Documentation: https://psap-cli.readthedocs.io.


1. Install psap
--------
.. code-block:: python

   python setup.py install


2. Annotate peptide sequences and save feature data frame to pickle.
--------
.. code-block:: python

   psap annotate -f /path/to/peptide.fasta -o /output/directory       

This step is required before using any dataset for training and prediction

3. Train RandomForest classifier and save to pickle.
--------
.. code-block:: python

   psap train -df /path/to/annotated/peptide/trainingset.pkl -o /output/directory  


4. Predict llps formation probability.
--------
.. code-block:: python

   psap predict -m /path/to/classifier.joblid -df /path/to/annotated/peptide/testset.pkl -o /output/directory


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
