# Define the truth values for angels and devils
ANGEL = True
DEVIL = False

# Define the statements
def mason_statement(mason, sofia):
    return mason == sofia

def sofia_statement(mason, sofia):
    return mason != sofia

# Check the possible scenarios
scenarios = [
    (ANGEL, ANGEL),
    (ANGEL, DEVIL),
    (DEVIL, ANGEL),
    (DEVIL, DEVIL)
]

# Find the valid scenario
for mason, sofia in scenarios:
    if (mason == ANGEL and mason_statement(mason, sofia)) or (mason == DEVIL and not mason_statement(mason, sofia)):
        if (sofia == ANGEL and sofia_statement(mason, sofia)) or (sofia == DEVIL and not sofia_statement(mason, sofia)):
            print(f"Mason is {'an angel' if mason else 'a devil'}, and Sofia is {'an angel' if sofia else 'a devil'}")