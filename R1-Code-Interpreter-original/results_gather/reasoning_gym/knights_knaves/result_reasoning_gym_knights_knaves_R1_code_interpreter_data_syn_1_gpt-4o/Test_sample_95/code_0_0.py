# Define the truth-telling nature of heros and villains
def is_hero(statement_truth):
    return statement_truth

def is_villain(statement_truth):
    return not statement_truth

# Victoria's statement: "Scarlett is a hero if and only if Scarlett is a villain."
victoria_statement = lambda scarlett_is_hero, scarlett_is_villain: scarlett_is_hero == scarlett_is_villain

# Scarlett's statement: "If Victoria is a hero then Victoria is a villain."
scarlett_statement = lambda victoria_is_hero, victoria_is_villain: not victoria_is_hero or victoria_is_villain

# Check all combinations
results = []
for victoria_is_hero in [True, False]:
    for scarlett_is_hero in [True, False]:
        victoria_is_villain = not victoria_is_hero
        scarlett_is_villain = not scarlett_is_hero
        
        victoria_truth = is_villain(victoria_statement(scarlett_is_hero, scarlett_is_villain))
        scarlett_truth = is_villain(scarlett_statement(victoria_is_hero, victoria_is_villain))
        
        if victoria_truth and scarlett_truth:
            results.append((victoria_is_hero, scarlett_is_hero))

print(results)