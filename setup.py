from pathlib import Path

from setuptools import find_packages
from setuptools import setup


requirements = ["numpy", "scipy", "matplotlib"]

current_directory = Path(__file__).parent

setup(
    author="Christopher Klugmann",
    author_email="ck@quality-match.com",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.9",
    ],
    description="A module to determine sample size or confidence intervals for proportions in multinomial populations.",
    install_requires=requirements,
    include_package_data=True,
    keywords="multinomial_estimation",
    name="multinomial_estimation",
    packages=find_packages(),
    url="https://github.com/cklugmann/multinomial-proportions",
    version="0.0.1",
    zip_safe=False,
)
