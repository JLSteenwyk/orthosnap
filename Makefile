install:
	# install
	python setup.py install

test.unit:
	# test units
	python -m pytest -m "not integration"

test.integration:
	# test end-to-end function
	python -m pytest --basetemp=output -m "integration"

test.fast:
	python -m pytest -m "not (integration or slow)"
	python -m pytest --basetemp=output -m "integration and not slow"

test:
	## execute both integration and unit tests
	# unit tests
	python -m pytest -m "not integration"
	# integration tests
	python -m pytest --basetemp=output -m "integration"


# used by GitHub actions during CI workflow
test.coverage: coverage.unit coverage.integration

coverage.unit:
	python -m pytest --cov=./ -m "not integration" --cov-report=xml:unit.coverage.xml

coverage.integration:
	rm tests/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.*
	python -m pytest --basetemp=output --cov=./ -m "integration" --cov-report=xml:integration.coverage.xml
	rm tests/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.*
