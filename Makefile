install:
	# install
	python setup.py install

test.unit:
	# test units
	python -m pytest -m "not integration"

test.integration:
	# test end-to-end function
	python -m pytest --basetemp=output -m "integration"

test:
	## execute both integration and unit tests
	# unit tests
	python -m pytest -m "not integration"
	# integration tests
	python -m pytest --basetemp=output -m "integration"

test.fast:
	python -m pytest -m "not (integration or slow)"
	python -m pytest --basetemp=output -m "integration and not slow"

# used by GitHub actions during CI workflow
test.coverage: coverage.unit coverage.integration

coverage.unit:
	python -m pytest --cov=./ -m "not integration" --cov-report=xml:unit.coverage.xml

coverage.integration:
	rm -rf output/
	mkdir output/
	python -m pytest --basetemp=output --cov=./ -m "integration" --cov-report=xml:integration.coverage.xml