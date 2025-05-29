def is_altruist(scarlett, daniel):
    # Scarlett's statement: If Scarlett is an altruist, then Daniel is an altruist
    scarlett_statement = not scarlett or daniel
    # Daniel's statement: Daniel is an altruist or Scarlett is an egoist
    daniel_statement = daniel or not scarlett
    return scarlett_statement and daniel_statement

# Check all possibilities
possibilities = [(True, True), (True, False), (False, True), (False, False)]
results = [(scarlett, daniel) for scarlett, daniel in possibilities if is_altruist(scarlett, daniel)]

print(results)