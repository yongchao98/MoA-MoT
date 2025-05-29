# Define the possible states
# True represents hero, False represents villain
def check_statements(benjamin, oliver):
    # Benjamin's statement: "Oliver is a villain or Benjamin is a hero"
    benjamin_statement = (not oliver) or benjamin
    
    # Oliver's statement: "If Oliver is a hero then Benjamin is a villain"
    oliver_statement = (not oliver) or (not benjamin)
    
    # Check if the statements are consistent with their nature
    benjamin_truth = benjamin == benjamin_statement
    oliver_truth = oliver == oliver_statement
    
    return benjamin_truth and oliver_truth

# Check all combinations
results = []
for benjamin in [True, False]:  # True for hero, False for villain
    for oliver in [True, False]:
        if check_statements(benjamin, oliver):
            results.append((benjamin, oliver))

print(results)