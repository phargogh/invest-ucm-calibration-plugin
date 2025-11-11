
env:
	# presumes python>=3.11 is available and on the PATH
	mamba create -p ./env $(python -c "import tomllib;print(' '.join(tomllib.load(open('pyproject.toml', 'rb'))['tool']['natcap']['invest']['conda_dependencies']))")
