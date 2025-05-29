# Define the truth-telling nature of heros and villains
def is_hero(statement_truth):
    return statement_truth

def is_villain(statement_truth):
    return not statement_truth

# Sofia's statement: "Isabella is a hero if and only if Isabella is a villain."
# This is a contradiction, so Sofia must be a villain.
sofia_statement_truth = False
sofia_is_villain = is_villain(sofia_statement_truth)

# Isabella's statement: "Sofia is a hero."
# Since Sofia is a villain, this statement is false, so Isabella must be a villain.
isabella_statement_truth = False
isabella_is_villain = is_villain(isabella_statement_truth)

print(f"Sofia is a {'villain' if sofia_is_villain else 'hero'}, and Isabella is a {'villain' if isabella_is_villain else 'hero'}.")