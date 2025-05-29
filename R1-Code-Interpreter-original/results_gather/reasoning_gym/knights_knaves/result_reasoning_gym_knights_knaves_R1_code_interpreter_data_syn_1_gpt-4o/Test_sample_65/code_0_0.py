# Define the truth values for sage and fool
SAGE = True
FOOL = False

# Define the statements
def mason_statement(mason, noah):
    return mason == SAGE or noah == SAGE

def noah_statement(mason, noah):
    return (mason == FOOL) == (noah == SAGE)

# Check all possible combinations
results = []
for mason in [SAGE, FOOL]:
    for noah in [SAGE, FOOL]:
        if mason_statement(mason, noah) == (mason == SAGE) and noah_statement(mason, noah) == (noah == SAGE):
            results.append((mason, noah))

# Print the results
print(results)