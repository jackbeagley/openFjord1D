import os
import re
import subprocess
import sys
from pathlib import Path

from skbuild import setup

from setuptools.command.develop import develop


# Convert distutils Windows platform specifiers to CMake -A arguments
PLAT_TO_CMAKE = {
    "win32": "Win32",
    "win-amd64": "x64",
    "win-arm32": "ARM",
    "win-arm64": "ARM64",
}

def initialise_submodules(dir = "."):
	subprocess.run(["git", "submodule", "update", "--init", "--recursive"], cwd=dir, check=True)

class CustomDevelopCommand(develop):
	def __init__(self, project_dir):
		develop.__init__(self)
		self.project_dir = project_dir

	def run(self):
		initialise_submodules(self.project_dir)
		develop.run(self)

package_dir = Path(__file__).resolve().parent
extern_dir = os.path.abspath(os.path.join(package_dir, 'extern'))

initialise_submodules(str(package_dir))

# The information here can also be placed in setup.cfg - better separation of
# logic and declaration, and simpler if you include description/version in a file.
setup(
    name="openFjord1D",
	packages=["openFjord1D"],
    version="1.0.0",
    author="Jack Beagley",
    author_email="jack.beagley@outlook.com",
#    description="A test project using pybind11 and CMake",
#    long_description="",
#    ext_modules=[CMakeExtension("openFjord1D")],
    cmake_source_dir="openFjord1D",
	cmake_args=[f"-DEXTERN_DIR={extern_dir}"],
#    cmdclass={
#		"develop": CustomDevelopCommand,
#		},
    zip_safe=False,
    python_requires=">=3.11",
	install_requires=[
		"numpy",
		],
)
