# Define the possible states
altruist = True
egoist = False

# Owen's statement: "Benjamin is an egoist and Owen is an altruist."
def owen_statement(owen, benjamin):
    return (benjamin == egoist) and (owen == altruist)

# Benjamin's statement: "Owen is an altruist if and only if Benjamin is an altruist."
def benjamin_statement(owen, benjamin):
    return (owen == altruist) == (benjamin == altruist)

# Check the scenarios
# Scenario 1: Owen is an altruist
owen_is_altruist = altruist
benjamin_is_altruist = egoist
scenario1 = owen_statement(owen_is_altruist, benjamin_is_altruist) and benjamin_statement(owen_is_altruist, benjamin_is_altruist)

# Scenario 2: Owen is an egoist
owen_is_altruist = egoist
benjamin_is_altruist = altruist
scenario2 = not owen_statement(owen_is_altruist, benjamin_is_altruist) and not benjamin_statement(owen_is_altruist, benjamin_is_altruist)

print("Scenario 1 (Owen is altruist):", scenario1)
print("Scenario 2 (Owen is egoist):", scenario2)