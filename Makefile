.PHONY: check, tcheck, pycodestyle, pyflakes, flake8, lint, wc, clean, clobber, upload

XARGS := xargs $(shell test $$(uname) = Linux && echo -r)

check:
	env PYTHONPATH=. python -m discover -v

tcheck:
	env PYTHONPATH=. trial --rterrors test

pytest:
	env PYTHONPATH=. pytest --verbose

pycodestyle:
	find bin civaclib test -name '*.py' -print0 | $(XARGS) -0 pycodestyle --ignore E402,W504

pyflakes:
	find bin civaclib test -name '*.py' -print0 | $(XARGS) -0 pyflakes

flake8:
	find bin civaclib test -name '*.py' -print0 | $(XARGS) -0 flake8 --ignore E402,W504

wc:
	find . -path './.tox' -prune -o -path './build' -prune -o -path './dist' -prune -o -name '*.py' -print0 | $(XARGS) -0 wc -l

clobber: clean
	rm -fr .tox
