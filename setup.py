from setuptools import setup, find_packages

def readme_file():
    with open("README.md", encoding="utf-8") as rf:
        return rf.read()

setup(
    name='MetaSAG',
    version='1.3.10',
    description='A microbial single-amplified genome droplet sequencing analysis toolkit',
    long_description=readme_file(),
    long_description_content_type="text/markdown",
    author='LiangCheng',
    author_email='liangcheng@hrbmu.edu.cn',
    url='https://github.com/liangcheng-hrbmu/MetaSAG',
    packages=find_packages(),
    package_data={
        'MetaSAG': [
            '*.R',
            '*.txt',
            '*.xlsx',
            'trimmomatic/*',
        ],
    },
    include_package_data=True,

    scripts=[
        'MetaSAG/k-mer.py',
        'MetaSAG/k_mer_union.py',
        'MetaSAG/MetaK_Lytic.py',
        'MetaSAG/phage_predictor.py',
        'MetaSAG/strain_analysis_multi_sgb.R'
    ],

    python_requires='>=3.8',

    install_requires=[
        'setuptools==59.8.0',
        'numpy',
        'pandas',
        'openpyxl',       
        'matplotlib',
        'seaborn',
        'scipy',          
        'umap-learn',     
        'scikit-learn',   
        'statsmodels',    
        'biopython',      
        'tensorflow>=2.0',
        'h5py',           
    ],
)