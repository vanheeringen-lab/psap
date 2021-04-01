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

CLI interface for the PSAP classifier. PSAP implements a RandomForest approach to predict the probability of proteins to mediate protein phase separation (PPS). Initially, a set of protein sequences is annotated with biochemical features wich are subsequently used to train a RandomForest (scikit-learn) classifier. The trained classifier is exported to json format and can be used to predict the llps class probability (PSAP_score) for new samples. The default model was trained on the human reference proteome `000005640_9606<ftp://ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/Eukaryota/UP000005640_9606.fasta.gz>`_ and is included in the installation.  

**Publication**
| Mierlo, G., Jansen, J. R. G., Wang, J., Poser, I., van Heeringen, S. J., & Vermeulen, M. (2021). Predicting protein condensate formation using machine learning. Cell Reports, 34(5), 108705. https://doi.org/10.1016/j.celrep.2021.108705.


* Free software: MIT license

========
Getting Started
========

1. *Install psap*
--------
.. code-block:: bash
   
   git clone https://github.com/vanheeringen-lab/psap.git
   cd psap && python setup.py install
   
2. *Train classifier*
--------
.. code-block:: python

   psap train -f /path/to/peptide-trainingset.fasta  -o /output/directory  
The trained RandomForest classifier is exported to json format and stored in the output directory.

3. *Predict llps score for peptide instances*
--------
.. code-block:: python

   psap predict -m /path/to/model.json -f /path/to/peptid-testset.fasta -o /output/directory
   
When no model is provided (-m) psap loads the default classifier stored in /data/model.

4. *Annotate petides (optional)*
--------
.. code-block:: python

   psap annotate -f /path/to/peptide.fasta -o /output/directory    

Annotates a peptide fasta with biochemical features. This step is included in train and predict.



Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
