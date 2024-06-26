"""
Script to run exPfact multiple times using a set of initial guesses (protection factors).
Filenames must be of the kind 'fname'+'{number}', with number = 1, ..., N, and formatted as specified by exPfact's documentation.
It is necessary to modify the script before running, specifying exPfact's input arguments.

Usage:
Set as working directory the directory in which data are (e.g., (..)/exPfact/data). From the terminal, run:

python ../python/run_multiple.py
"""

import subprocess

N = 100

guess = "--pfact"

for n in range(1, N+1):
	input_file = f"guess_{n}.pfact"
	subprocess.run(['python',
					'../python/exPfact.py',
					"--temp", "277.15",
					"--pH", "7",
					"--dexp", "gb1.dexp",
					"--ass", "gb1.list",
					"--seq", "gb1.seq",
					"--out", f"out_{n}",
					guess, f"guess_{n}.pfact"
					])
