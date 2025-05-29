# Define the possible states
sages = {"Logan": False, "Jacob": False}

# Check Case 1: Jacob is a sage
sages["Jacob"] = True
sages["Logan"] = True  # Because Jacob says Logan is a sage

# Verify Logan's statement: If Jacob is a sage, then Logan is a sage
logans_statement = not sages["Jacob"] or sages["Logan"]

# Check if both statements hold true
result = logans_statement and sages["Jacob"] == True

print(result)