# Publishing SuppOrth

The version is defined in `src/supporth/__init__.py`. Use a new version for every
upload; PyPI does not allow a published distribution to be replaced.

## Build and check

Use an existing Python interpreter in a virtual environment:

```sh
python3.10 -m venv /tmp/supporth-release
source /tmp/supporth-release/bin/activate
python -m pip install '.[dev,release]'
python -m pytest tests
python -m build
python -m twine check --strict dist/*
```

Build from a clean checkout so `dist/` contains only the intended release.
`python -m build` creates both the source archive and a wheel built from that
archive. Check installation in another environment, outside the checkout:

```sh
python3.10 -m venv /tmp/supporth-wheel-check
/tmp/supporth-wheel-check/bin/python -m pip install dist/supporth-0.4.1-py3-none-any.whl
/tmp/supporth-wheel-check/bin/python -m pip check
/tmp/supporth-wheel-check/bin/suppOrth --version
/tmp/supporth-wheel-check/bin/supporth run --help
```

These checks validate the Python distribution. Predictor installation and live
orthology runs require separate checks with the upstream software and Pfam data.

## Upload

Create PyPI and TestPyPI accounts and use their API tokens when Twine asks for
credentials. Tokens must not be stored in the repository.

```sh
python -m twine upload --repository testpypi dist/*
python -m twine upload dist/*
```

TestPyPI and PyPI are separate services. Confirm the TestPyPI metadata and files
before uploading to PyPI. Stable releases are installed with
`python -m pip install supporth`.

The package name must be available or owned by your PyPI account. An apparent
absence from the public index does not guarantee that PyPI will allow the name.
