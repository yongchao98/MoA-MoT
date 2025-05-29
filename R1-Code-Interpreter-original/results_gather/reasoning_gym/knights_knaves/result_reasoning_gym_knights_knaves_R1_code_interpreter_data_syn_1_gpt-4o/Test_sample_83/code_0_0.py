# Define the possible states
altruist = True
egoist = False

# Zoey's statement: If Ella is an altruist, then Zoey is an altruist
def zoey_statement(ella, zoey):
    return not ella or zoey

# Ella's statement: Zoey is an altruist if and only if Zoey is an egoist
def ella_statement(zoey):
    return zoey == (not zoey)

# Check all combinations
results = []
for zoey in [altruist, egoist]:
    for ella in [altruist, egoist]:
        if zoey_statement(ella, zoey) and not ella_statement(zoey):
            results.append((zoey, ella))

# Print the result
print(results)