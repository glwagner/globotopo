from setuptools import setup

setup(
    name = 'globotopo',
    version = '0.1',
    description = "A tool for loading and cutting topographic data.", 
    url = 'http://github.com/glwagner/globotopo',
    author = 'Gregory L. Wagner',
    author_email = 'wagner.greg@gmail.com',
    license = 'MIT',
    packages = ['globotopo'],
    install_requires = [
        'numpy', 
        'matplotlib', 
    ],
    zip_safe = False,
)
