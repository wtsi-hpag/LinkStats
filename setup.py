import sys
from subprocess import Popen

import numpy as np
from setuptools import find_packages, setup
from setuptools.command.install import install as setup_install

VERSION = "0.0.1"

LINKSTATS_C = "_LinkStats_C"


def install_linkstats_c():
    for name, cmd in (
        (
            "setting meson version",
            f"meson rewrite kwargs set project / version {VERSION}",
        ),
        (
            "meson setup",
            f"meson setup --buildtype=release --unity on --prefix={sys.prefix} builddir",
        ),
        ("_LinkStats_C compile", "meson compile -C builddir"),
        ("_LinkStats_C test", "meson test -C builddir"),
        ("_LinkStats_C install", "meson install -C builddir"),
    ):
        with Popen(
            cmd.split(),
            cwd=LINKSTATS_C,
            stdout=sys.stdout,
            stderr=sys.stderr,
        ) as proc:
            if proc.wait() != 0:
                sys.exit("Error during " + name)


class LinkStatsInstall(setup_install):
    def run(self):
        install_linkstats_c()
        setup_install.run(self)


def main():

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
            "tqdm>=4.61.1",
            "networkx>=2.6.3",
        ),
        description="Collect stats from aligned linked-reads",
        entry_points="""
              [console_scripts]
              LinkStats=LinkStats.main:cli
          """,
        cmdclass={"install": LinkStatsInstall},
    )


if __name__ == "__main__":
    main()
