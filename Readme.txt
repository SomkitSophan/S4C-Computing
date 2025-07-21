Publisher: CSSRG-LAB of the King Mongkut's Institute of Technology Ladkrabang (KMITL), Thailand.

Requirements: Python3.x, GPS data of RINEX2 (for example file in "positioning/Input/").

Simulation steps:
       1. Copy a folder of "positioning" to your computer at the path of "C:/" (Windows) or "~/home/" (Ubuntu).
       2. Install important packges:
		 => PyCharm (Community edition):
			1) "numpy"
			2) "pandas"
			3) "matplotlib"
			4) "pyunpack"
			5) "patool"
		 => Ubuntu 20.04 on Win10 ("root" user):
			1) "sudo apt update"
			2) "sudo apt install python3-pip"
			3) "sudo apt-get install python3-tk"
			4) "pip install pandas"
			5) "pip install matplotlib"
			6) "pip install pyunpack"
			7) "pip install patool"
			8) "pip install paramiko"
	3. Define "Input_path" & "Output_path" along with your computer (line 21,22), and a station code of RINEX file (line 27).
	4. Copy your RINEX files, "*.yyo" files to "~/positioning/Input/obs/" and "*.yyn" to "~/positioning/Input/nav/".
	5. Press button "Run" (PyCharm) or enter "/usr/bin/python3 /home/positioning/*.py" (Ubuntu, root user)

Simulation results:
	1. Each error is computed by the sampling time of 30 seconds. Then the 10 results (5 minutes) are averaged and logged as the "*.txt" file into the directory of "~/positioning/Output/".
	2. Definitions of each parameter in the "*.txt" file are explained in the header for the daily file (PositioningErr_daily_V1.py), no a near-real time file (PositioningErr_rt_V1.py). For example file is in "positioning/Output/"