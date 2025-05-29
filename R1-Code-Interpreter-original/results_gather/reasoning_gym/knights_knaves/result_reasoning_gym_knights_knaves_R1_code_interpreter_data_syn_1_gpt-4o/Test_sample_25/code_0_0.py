# Define the possible states
sages = {"Lucas": False, "Jacob": False}  # Assume both are fools initially

# Check Lucas's statement: "Lucas is a sage if and only if Jacob is a sage"
lucas_statement = (sages["Lucas"] == sages["Jacob"])

# Check Jacob's statement: "If Jacob is a sage then Lucas is a fool"
jacob_statement = (not sages["Jacob"] or not sages["Lucas"])

# Both statements should be true for the scenario to be consistent
consistent_scenario = lucas_statement and jacob_statement

print(consistent_scenario)