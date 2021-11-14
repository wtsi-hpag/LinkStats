import os
import sys
from glob import glob
from pathlib import Path
from subprocess import Popen

import numpy as np
from setuptools import Extension, find_packages, setup

VERSION = "0.0.1"
TMP_PREFIX = Path(os.getcwd()).absolute() / "libLinkStats_Prefix"


def build_liblinkstats():
    with Popen(
        f"meson rewrite kwargs set project / version {VERSION}".split(),
        cwd="libLinkStats",
        stdout=sys.stdout,
        stderr=sys.stderr,
    ) as proc:
        if proc.wait() != 0:
            sys.exit("Error setting meson version")

    with Popen(
        f"meson setup --buildtype=release --unity on --prefix={TMP_PREFIX} builddir".split(),
        cwd="libLinkStats",
        stdout=sys.stdout,
        stderr=sys.stderr,
    ) as proc:
        if proc.wait() != 0:
            sys.exit("Error during meson setup")

    with Popen(
        f"meson install -C builddir".split(),
        cwd="libLinkStats",
        stdout=sys.stdout,
        stderr=sys.stderr,
    ) as proc:
        if proc.wait() != 0:
            sys.exit("Error during libLinkStats compile")

    return tuple(glob(str(TMP_PREFIX / "**" / "libLinkStats.*"), recursive=True))


def main():

    libLinkStats_object = build_liblinkstats()
    assert len(libLinkStats_object) == 1

    setup(
        name="LinkStats",
        version=VERSION,
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
        ext_modules=[
            Extension(
                "_LinkStats_C",
                sources=["LinkStats_PyCExt.cpp"],
                include_dirs=[TMP_PREFIX / "include", np.get_include()],
                library_dirs=[os.getcwd()],
                extra_compile_args=["-Ofast"],
                extra_objects=libLinkStats_object,
                extra_link_args=["-lstdc++"],
                language="c++20",
            )
        ],
        entry_points="""
              [console_scripts]
              LinkStats=LinkStats.main:cli
          """,
    )


if __name__ == "__main__":
    main()
