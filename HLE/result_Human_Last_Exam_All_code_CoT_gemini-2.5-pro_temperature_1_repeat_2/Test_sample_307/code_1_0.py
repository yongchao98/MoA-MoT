def predict_drosophila_outcome(mother_cholesterol_mg_per_L, offspring_cholesterol_mg_per_L):
    """
    Analyzes Drosophila development based on dietary rules.

    Args:
        mother_cholesterol_mg_per_L: The cholesterol concentration in the mother's diet.
        offspring_cholesterol_mg_per_L: The cholesterol concentration in the offspring's diet.

    Returns:
        A string describing the developmental outcome for the offspring.
    """
    
    # Rule: A diet is insufficient if cholesterol is less than 250mg/L.
    # The problem states larvae cannot survive to adulthood on 2mg/L or 0mg/L cholesterol.
    # We generalize this to any amount significantly less than the "normal" 250mg/L.
    if offspring_cholesterol_mg_per_L < 250:
        # This is the most specific outcome for developmental failure due to lack of molting hormone.
        outcome = "No eclosion to adulthood"
    else:
        outcome = "Normal development"

    # Note on the mother's diet:
    # A mother reared on 0mg/L cholesterol (i.e., pure cholestanol) would not survive to reproduce.
    # The code proceeds assuming the question wants an analysis of the offspring's diet regardless.
    if mother_cholesterol_mg_per_L < 250:
        mother_viability = "non-viable to be reared on"
    else:
        mother_viability = "viable"
        
    return outcome, mother_viability

# --- Scenario 1 Analysis ---
# Mothers' diet: 250mg/L cholesterol.
# Offspring's diet: 250mg/L cholestanol (which means 0mg/L cholesterol).
scenario1_mother_diet = 250
scenario1_offspring_diet = 0
outcome1, _ = predict_drosophila_outcome(scenario1_mother_diet, scenario1_offspring_diet)


# --- Scenario 2 Analysis ---
# Mothers' diet: 250mg/L cholestanol (0mg/L cholesterol).
# Offspring's diet: 250mg/L cholestanol and 2mg/L cholesterol.
scenario2_mother_diet = 0
scenario2_offspring_diet = 2
outcome2, mother_viability_s2 = predict_drosophila_outcome(scenario2_mother_diet, scenario2_offspring_diet)

# --- Print the results ---
print("--- Analysis Report ---")
print("Question 1: What will happen to Drosophila fed on a diet of 250mg/L cholestanol when raised from eggs from mothers that were reared on 250mg/L cholesterol?")
print(f"Offspring diet contains {scenario1_offspring_diet}mg/L of usable cholesterol.")
print(f"Predicted Outcome 1: {outcome1}\n")

print("Question 2: What will happen to Drosophila fed on a diet of 250mg/L cholestanol and 2mg/L cholesterol when raised from eggs from mothers that were reared on 250mg/L cholestanol?")
print(f"Note: A mother's diet of {scenario2_mother_diet}mg/L cholesterol is {mother_viability_s2}.")
print(f"Offspring diet contains {scenario2_offspring_diet}mg/L of usable cholesterol, which is insufficient.")
print(f"Predicted Outcome 2: {outcome2}\n")

print("--- Final Answer ---")
print(f"The combined answer is: ({outcome1}, {outcome2})")