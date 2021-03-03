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