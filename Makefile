.PHONY: fmt
fmt:
	black scripts tests
	flake8 scripts tests

.PHONY: unit-test
unit-test:
	python3 -m unittest discover --start-directory tests/unit --top-level-directory . --verbose

help: ## display help for this makefile
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 
'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
