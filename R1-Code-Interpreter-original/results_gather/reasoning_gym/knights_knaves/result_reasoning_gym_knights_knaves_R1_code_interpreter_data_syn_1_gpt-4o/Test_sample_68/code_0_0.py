# Define the truth-telling nature of angels and devils
def is_angel(statement_truth):
    return statement_truth

def is_devil(statement_truth):
    return not statement_truth

# Ethan's statement: "If Oliver is a devil then Oliver is an angel"
def ethan_statement(oliver_is_devil):
    return oliver_is_devil and not oliver_is_devil

# Oliver's statement: "Ethan is an angel if and only if Oliver is an angel"
def oliver_statement(ethan_is_angel, oliver_is_angel):
    return ethan_is_angel == oliver_is_angel

# Check all possibilities
results = []
for ethan_is_angel in [True, False]:
    for oliver_is_angel in [True, False]:
        ethan_is_devil = not ethan_is_angel
        oliver_is_devil = not oliver_is_angel
        
        # Evaluate statements
        ethan_truth = ethan_statement(oliver_is_devil)
        oliver_truth = oliver_statement(ethan_is_angel, oliver_is_angel)
        
        # Check consistency with their nature
        if is_angel(ethan_is_angel) == ethan_truth and is_angel(oliver_is_angel) == oliver_truth:
            results.append((ethan_is_angel, oliver_is_angel))

print(results)