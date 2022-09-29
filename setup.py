import setuptools
from dna_seq_analysis.__init__ import __VERSION__, __STEPS__

with open("README.md", "r") as fh:
    long_description = fh.read()

with open('requirements.txt') as fp:
    install_requires = fp.read()

entrys = ['dna-seq=dna_seq_analysis.dna_seq_analysis:main',]
for assay in __STEPS__:
    entrys.append(f'dna-seq_{assay}=dna_seq_analysis.{assay}.{assay}:main')
entrys.append(f'dna-seq_multi=dna_seq_analysis.multi.multi:main')
entry_dict = {
        'console_scripts': entrys,
}


setuptools.setup(
    name="dna_seq_analysis",
    version=__VERSION__,
    author="wuqi",
    author_email="wuqi@singleronbio.com",
    description="Dna Seq Analysis Pipelines",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/singleron-RD/Dna-seq_Analysis",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    include_package_data=True,
    entry_points=entry_dict,
    install_requires=install_requires,
)