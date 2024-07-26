# Example code for running nf-test. It should be updated to include other executors

## Run tests in tests (only for workflow)
nf-test test --profile docker 

## Run tests in modules
nf-test test modules/ --profile docker 