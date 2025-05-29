# Define the truth-telling nature of heroes and villains
def is_hero(statement_truth):
    return statement_truth

def is_villain(statement_truth):
    return not statement_truth

# Abigail's statement: "Mason is a villain and Mason is a hero."
abigail_statement_truth = False  # This statement is a contradiction, hence false

# Mason's statement: "Abigail is a hero."
mason_statement_truth = False  # Since Abigail's statement is false, Mason's statement is also false

# Determine the nature of Abigail and Mason
abigail_is_hero = is_hero(abigail_statement_truth)
mason_is_hero = is_hero(mason_statement_truth)

# Print the results
print(f"Abigail is a {'hero' if abigail_is_hero else 'villain'}, and Mason is a {'hero' if mason_is_hero else 'villain'}.")