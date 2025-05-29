# Define the possible states
# True represents a pioneer, and False represents a laggard
def check_statements(victoria, amelia):
    # Victoria's statement: If Amelia is a pioneer, then Amelia is a laggard
    victoria_statement = not amelia or not amelia
    
    # Amelia's statement: Victoria is a pioneer if and only if Amelia is a pioneer
    amelia_statement = (victoria == amelia)
    
    # Both statements must be consistent with their roles
    return (victoria == victoria_statement) and (amelia == amelia_statement)

# Check all possible combinations
results = []
for victoria in [True, False]:
    for amelia in [True, False]:
        if check_statements(victoria, amelia):
            results.append((victoria, amelia))

# Print the results
print(results)