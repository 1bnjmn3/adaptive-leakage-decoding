# Local tasks for the heralded-leakage decoding study.
# Recipe lines start with a tab. GNU make is required.
# Run `make venv` once, then `make test`, `make bench-quick`, or `make hw`.

.PHONY: venv test bench bench-quick hw zip clean

venv:
	python3 -m venv .venv && .venv/bin/pip install -e ".[test]"

test:
	.venv/bin/python -m pytest -q

bench:
	.venv/bin/leakdecode-bench

bench-quick:
	.venv/bin/leakdecode-bench --quick

hw:
	make -C hw

zip:
	cd .. && zip -r riverlane-project.zip riverlane -x 'riverlane/.venv/*' 'riverlane/**/__pycache__/*' 'riverlane/hw/build/*' 'riverlane/*.egg-info/*' 'riverlane/src/*.egg-info/*'

clean:
	rm -rf .venv .pytest_cache hw/build results-ci src/*.egg-info *.egg-info
	find . -depth -type d -name __pycache__ -exec rm -rf {} +
	rm -f ../riverlane-project.zip
