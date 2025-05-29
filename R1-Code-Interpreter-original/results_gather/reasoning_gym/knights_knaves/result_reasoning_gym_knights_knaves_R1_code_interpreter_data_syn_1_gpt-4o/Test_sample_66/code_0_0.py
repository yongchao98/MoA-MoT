def check_statements(lily_is_hero, jack_is_hero):
    # Lily's statement: "Lily is a hero or Jack is a villain."
    lily_statement = lily_is_hero or not jack_is_hero
    
    # Jack's statement: "If Jack is a hero then Lily is a villain."
    jack_statement = not jack_is_hero or not lily_is_hero
    
    # Check if the statements are consistent with their roles
    lily_consistent = lily_statement if lily_is_hero else not lily_statement
    jack_consistent = jack_statement if jack_is_hero else not jack_statement
    
    return lily_consistent and jack_consistent

# Check both scenarios
scenario_1 = check_statements(lily_is_hero=True, jack_is_hero=False)
scenario_2 = check_statements(lily_is_hero=False, jack_is_hero=True)

print("Scenario 1 (Lily is a hero, Jack is a villain):", scenario_1)
print("Scenario 2 (Lily is a villain, Jack is a hero):", scenario_2)