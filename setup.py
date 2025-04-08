from setuptools import setup, find_packages

setup(
    name='ccm_demo',
    version='0.1.0',
    description="Python package for demonstrating simple proteomics integrations and one liners",
    author='Alper Celik',
    author_email='alper.celik@sickkids.ca',
    packages=find_packages(),
    zip_safe=False,
    package_data={"": ["*.json"]},
    scripts=["ccm_demo/scripts/run_mmseqs.sh"],
    include_package_data=True
)
