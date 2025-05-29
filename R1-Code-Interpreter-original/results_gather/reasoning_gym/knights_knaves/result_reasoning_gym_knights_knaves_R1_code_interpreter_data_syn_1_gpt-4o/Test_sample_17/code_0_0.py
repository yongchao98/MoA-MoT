def check_statements(sophia_is_altruist, oliver_is_altruist):
    # Sophia's statement: "Sophia is an altruist if and only if Oliver is an egoist."
    sophia_statement = (sophia_is_altruist == (not oliver_is_altruist))
    
    # Oliver's statement: "Sophia is an egoist if and only if Oliver is an altruist."
    oliver_statement = ((not sophia_is_altruist) == oliver_is_altruist)
    
    return sophia_statement and oliver_statement

# Check both scenarios
scenario_1 = check_statements(True, False)  # Sophia is altruist, Oliver is egoist
scenario_2 = check_statements(False, True)  # Sophia is egoist, Oliver is altruist

print("Scenario 1 (Sophia is altruist, Oliver is egoist):", scenario_1)
print("Scenario 2 (Sophia is egoist, Oliver is altruist):", scenario_2)