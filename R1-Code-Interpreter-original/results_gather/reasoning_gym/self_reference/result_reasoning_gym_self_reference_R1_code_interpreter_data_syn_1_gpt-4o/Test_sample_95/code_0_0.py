# Possible scenarios
possible_scenarios = 0

# Check for 5 true statements
if 5 in [2, 3, 5, 7]:  # Prime numbers
    if 2 in [4, 6, 8]:  # Composite numbers
        possible_scenarios += 1

# Check for 3 true statements
if 3 in [2, 3, 5, 7]:  # Prime numbers
    if 4 in [4, 6, 8]:  # Composite numbers
        possible_scenarios += 1

print(possible_scenarios)