from setuptools import find_packages, setup


def main():
    setup(
        name="LinkStats",
        version="0.0.1",
        packages=find_packages(),
        include_package_data=True,
        author="Ed Harry",
        author_email="edward.harry@sanger.ac.uk",
        python_requires=">=3.8.10",
        install_requires=(
            "click>=8.0.1",
            "pandas>=1.2.5",
            "numpy>=1.20.3",
            "seaborn>=0.11.1",
            "matplotlib>=3.4.2",
            "pysam>=0.16.0.1",
            "tqdm>=4.61.1",
            "networkx>=2.6.3",
        ),
        description="Collect stats from aligned linked-reads",
        entry_points="""
              [console_scripts]
              LinkStats=LinkStats.main:cli
          """,
    )


if __name__ == "__main__":
    main()
