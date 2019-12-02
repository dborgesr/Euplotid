from distutils.core import setup
setup(
  name = 'euplotid',
  packages = ['euplotid'], 
  version = '0.2',
  description = 'A linux-based platform to identify, predict, and assess the difficulty of inducing transition mutations on Cis-Regulatory Elements constrained within Insulated Neighborhoods',
  author = 'Diego Borges',
  author_email = 'dborgesr@mit.edu',
  url = 'https://github.com/dborgesr/Euplotid', 
  download_url = 'https://github.com/dborgesr/Euplotid/archive/0.2.tar.gz', 
  keywords = ['CRISPR','neural_network', 'docker', 'linux'], 
  install_requires=[
    'readline','jupyter','plotly','sense-hat','mygene','myvariant','pysam','networkx','pandas','joblib','biopython','pyliftover','scipy','louvain','nbpresent',
  ],
  classifiers = [],
)
