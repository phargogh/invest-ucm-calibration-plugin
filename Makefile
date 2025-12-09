
.PHONY: coverage,env,test

test:
	python -m pytest -n auto -vs tests/tests.py --pdb

coverage:
	coverage run --source=src,env/lib/python3.13/site-packages/invest_ucm_calibration -m pytest tests/tests.py
	coverage html

env:
	# presumes python>=3.11 is available and on the PATH
	mamba create -p ./env $(python -c "import tomllib;print(' '.join(tomllib.load(open('pyproject.toml', 'rb'))['tool']['natcap']['invest']['conda_dependencies']))")

