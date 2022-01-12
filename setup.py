# Copyright (c) 2022 Ed Harry, Wellcome Sanger Institute, Genome Research Limited
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import sys
from pathlib import Path
from shutil import copyfile
from subprocess import PIPE, Popen

from setuptools import find_packages, setup
from setuptools.command.install import install as setup_install

with Popen("git describe".split(), stdout=PIPE, stderr=sys.stderr, text=True) as proc:
    if proc.wait() != 0:
        try:
            with open("version", "r") as vf:
                VERSION = next(vf).strip().replace("-", "_")
        except:
            sys.exit(
                "Could not determine version from 'git describe' or 'version' file"
            )
    else:
        VERSION = next(proc.stdout).strip().replace("-", "_")


LINKSTATS_C = "_LinkStats_C"


def install_linkstats_c():
    copyfile(Path(LINKSTATS_C) / "meson.build.in", Path(LINKSTATS_C) / "meson.build")

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
            cmd.split(), cwd=LINKSTATS_C, stdout=sys.stdout, stderr=sys.stderr,
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
        python_requires=">=3.8.0",
        install_requires=(
            "click>=8.0.1",
            "pandas>=1.2.5",
            "numpy>=1.20.3",
            "seaborn>=0.11.1",
            "matplotlib>=3.4.2",
            "tqdm>=4.61.1",
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
