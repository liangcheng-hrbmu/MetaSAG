from setuptools import setup, find_packages
#from distutils.core import setup

def readme_file():
    with open("README.md", encoding="utf-8") as rf:
        return rf.read()

setup(
    name='MetaSAG',
    version='1.0.20',
    description='Just Test PyPI Package',
    long_description=readme_file(),
    long_description_content_type="text/markdown",
    author='DMY',
    author_email='2023020560@hrbmu.edu.cn',
    url='https://github.com',
    packages=find_packages(),
    package_data={'':['READMESource/*']},
    include_package_data=True,
    #py_modules=['MetaSAG.BarcodeDeal','MetaSAG.BCFilter','MetaSAG.BinQCAnno','MetaSAG.CellHGT','MetaSAG.HUMAnNPath','MetaSAG.MetaPhlAnNAsign','SNPStrain'],
    python_requires='>=3.6'
)

