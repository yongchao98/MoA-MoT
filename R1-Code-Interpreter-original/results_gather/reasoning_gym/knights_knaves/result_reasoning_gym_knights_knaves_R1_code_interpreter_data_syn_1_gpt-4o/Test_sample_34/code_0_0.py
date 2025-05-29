def check_statements(luke_is_angel, penelope_is_angel):
    # Luke's statement: "Penelope is a devil"
    luke_statement = not penelope_is_angel
    
    # Penelope's statement: "If Penelope is an angel then Luke is a devil"
    penelope_statement = not penelope_is_angel or not luke_is_angel
    
    # Check if both statements are consistent with their nature
    luke_truth = luke_is_angel == luke_statement
    penelope_truth = penelope_is_angel == penelope_statement
    
    return luke_truth and penelope_truth

# Check both scenarios
scenario_1 = check_statements(True, False)  # Luke is an angel, Penelope is a devil
scenario_2 = check_statements(False, True)  # Luke is a devil, Penelope is an angel

print("Scenario 1 (Luke is an angel, Penelope is a devil):", scenario_1)
print("Scenario 2 (Luke is a devil, Penelope is an angel):", scenario_2)