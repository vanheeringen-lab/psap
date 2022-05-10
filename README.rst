========
psap
========


.. image:: https://github.com/vanheeringen-lab/psap/actions/workflows/python-app.yml/badge.svg
   :target:  https://github.com/vanheeringen-lab/psap

.. image:: https://github.com/vanheeringen-lab/psap/actions/workflows/continuous-deployment.yml/badge.svg
   :target:  https://github.com/vanheeringen-lab/psap

.. image:: https://badge.fury.io/py/psap.svg
   :target:  https://pypi.org/project/psap/

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
   :target:  https://anaconda.org/bioconda/psap
 
CLI interface for the PSAP classifier. PSAP implements a RandomForest approach to predict the probability of proteins to mediate protein phase separation (PPS). Initially, a set of protein sequences is annotated with biochemical features which are subsequently used to train a RandomForest (scikit-learn) classifier. The trained classifier is afterwards exported to json format and can be used to predict the llps class probability (PSAP_score) for unseen protein sequences. 

The default model was trained on the `human reference proteome <ftp://ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/Eukaryota/UP000005640_9606.fasta.gz>`_ with a list of literature curated PPS proteins for positive class labeling. Both can be found in /data.   

**Publication**
| Mierlo, G., Jansen, J. R. G., Wang, J., Poser, I., van Heeringen, S. J., & Vermeulen, M. (2021). Predicting protein condensate formation using machine learning. Cell Reports, 34(5), 108705. https://doi.org/10.1016/j.celrep.2021.108705.


* Free software: MIT license

================
Getting Started
================

1. *Install psap*
----------------------
**pip**

.. code-block:: bash
   
   pip install psap

**conda**

.. code-block:: bash
   
   conda install bioconda::psap
   
2. *Train classifier*
-----------------------
.. code-block:: python

   psap train -f /path/to/peptide-trainingset.fasta -l /path/to/known/pps-proteins.txt (optional) -o /path/to/output/directory (optional)
      
The trained RandomForest classifier will be exported to json format and stored in the output directory.

3. *Predict llps score for new peptide instances*
-----------------------------------------------
.. code-block:: python

   psap predict -f /path/to/peptide-testset.fasta -m /path/to/model.json (optional) -o /path/to/output/directory (optional)
   
psap loads the default classifier stored in `/data/model` when no model is provided with `-m`.
An example peptide fasta file can be found in the `psap/data/testset` folder. 

4. *Annotate peptide sequences (optional)*
---------------------------------
.. code-block:: python

   psap annotate -f /path/to/peptide.fasta -l /path/to/known/pps-proteins.txt (optional) -o /path/to/output/directory (optional)    

Manually annotate a peptide fasta with biochemical features. This step is included in train and predict.



Credits
-------

This package was adapted from the cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
