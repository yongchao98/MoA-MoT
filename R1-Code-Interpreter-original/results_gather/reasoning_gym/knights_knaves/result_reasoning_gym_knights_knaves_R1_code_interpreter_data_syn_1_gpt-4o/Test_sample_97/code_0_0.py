# Define the truth values for pioneers and laggards
pioneer = True
laggard = False

# Logan's statement: "Joseph is a pioneer and Logan is a pioneer."
def logan_statement(joseph, logan):
    return joseph == pioneer and logan == pioneer

# Joseph's statement: "Logan is a pioneer if and only if Logan is a laggard."
def joseph_statement(logan):
    return logan == pioneer if logan == laggard else False

# Check all possible combinations
results = []
for joseph in [pioneer, laggard]:
    for logan in [pioneer, laggard]:
        if logan_statement(joseph, logan) == (logan == pioneer) and joseph_statement(logan) == (joseph == pioneer):
            results.append((logan, joseph))

print(results)