#
# TESTING
# Integration testing
#   

BIN_DIR ?= $(abspath $(CURDIR)/../bin)
SELFTEST_DIR = $(wildcard $(CURDIR)/self-test)
EXAMPLE_SCRIPTS = $(sort $(wildcard $(CURDIR)/examples/**/test.sh) $(wildcard $(CURDIR)/examples/**/**/test.sh))

.PHONY: all
all: TEST_HEADER $(EXAMPLE_SCRIPTS) $(CONVERTCFG_SCRIPTS) UNIT_TEST 

.PHONY: clean
clean: $(sort $(wildcard $(CURDIR)/examples/**/out))
	@ $(RM) -r $^

.PHONY: TEST_HEADER
TEST_HEADER:
	@echo INTEGRATION TESTING:

.PHONY: $(EXAMPLE_SCRIPTS) $(CONVERTCFG_SCRIPTS)
$(EXAMPLE_SCRIPTS):
	@echo Running test $(notdir $(realpath $(dir $@)))
	@cd $(dir $@); \
	 mkdir -p out; \
	 $@ 1>out/stdout.log

.PHONY: UNIT_TEST
UNIT_TEST:
	@echo Running self test
	@cd $(SELFTEST_DIR); $(BIN_DIR)/mlp self-test
	@cd $(SELFTEST_DIR); $(BIN_DIR)/mlp self-test-dev

