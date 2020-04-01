from distutils.core import setup
setup(
  name = 'euplotid',
  packages = ['euplotid'], 
  version = '0.2',
  description = 'A quantized geometric model of the eukaryotic cell',
  author = 'Diego Borges',
  author_email = 'dborgesr@mit.edu',
  url = 'https://github.com/dborgesr/Euplotid', 
  download_url = 'https://github.com/dborgesr/Euplotid/archive/0.2.tar.gz', 
  keywords = ['CRISPR','neural_network', 'docker', 'linux'], 
  install_requires=[
    'readline','jupyter','plotly','sense-hat','mygene','myvariant','pysam','networkx','pandas','joblib','biopython','pyliftover','scipy','louvain',
  ],
  classifiers = [],
)
