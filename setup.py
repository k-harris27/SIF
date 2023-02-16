import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='SIF',
    version='0.0.2',
    author='Kieran Harris',
    author_email='kieran.harris@postgrad.manchester.ac.uk',
    description='A package containing a number of tools useful for modifying simulation input files.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/k-harris27/SIF',
    project_urls = {
    },
    license='MIT',
    packages=['SIF'],
    install_requires=[],
)
