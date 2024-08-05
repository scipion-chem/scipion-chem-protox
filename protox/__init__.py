# **************************************************************************
# *
# * Authors: Verónica Gamo (veronica.gamoparejo@usp.ceu.es)
# *
# * Biocomputing Unit, CNB-CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

# Scipion em imports
import os, subprocess
from subprocess import run
from scipion.install.funcs import InstallHelper

# Scipion chem imports
import pwchem

# Plugin imports
from .constants import PLUGIN_VERSION

_version_ = PLUGIN_VERSION
_logo = ""
_references = ['']

class Plugin(pwchem.Plugin):
	@classmethod
	def _defineVariables(cls):
		""" Return and write a variable in the config file. """

	@classmethod
	def defineBinaries(cls, env):
		""" Install the necessary packages. """

	########################### PACKAGE FUNCTIONS ###########################

	@classmethod
	def runProtox(cls, smiles_list, model, output_file, cwd=None):
			""" Run Protox command for given SMILES and model. """
			path= os.path.join(cls.getScriptsPath(), "protox3_api.py")
			program = f" python3 {path}"
			smiles_combined = ', '.join(smiles_list)
			args = f"-t smiles -m \"{model}\" -o {output_file} \"{smiles_combined}\""
			full_program = f'{program} {args}'
			run(full_program, env=cls.getEnviron(), shell=True)

	@classmethod
	def getPluginH(cls, path=""):
		import protox
		fnDir = os.path.split(protox.__file__)[0]
		return os.path.join(fnDir, path)
	
	@classmethod
	def getScriptsPath(cls):
		return cls.getPluginH('scripts/')
