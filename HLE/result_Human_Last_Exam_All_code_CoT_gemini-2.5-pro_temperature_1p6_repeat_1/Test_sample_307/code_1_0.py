def predict_drosophila_outcome(maternal_diet, offspring_diet):
    """
    Analyzes Drosophila development based on dietary rules.

    Args:
        maternal_diet (dict): A dictionary describing the mother's diet.
        offspring_diet (dict): A dictionary describing the offspring's diet.

    Returns:
        str: The predicted developmental outcome.
    """
    # Survival threshold for cholesterol is 250mg/L.
    survival_threshold_cholesterol = 250
    
    # Get the effective cholesterol from the offspring's diet. Cholestanol doesn't count.
    offspring_cholesterol = offspring_diet.get('cholesterol', 0)

    # According to the text, any diet less than 250mg/L cholesterol is insufficient
    # for larvae to survive to adulthood (e.g., 2mg/L or 0mg/L).
    # The ultimate fate is determined by the larval diet. Maternal provisioning
    # only provides a temporary buffer.
    if offspring_cholesterol >= survival_threshold_cholesterol:
        return "Normal development"
    else:
        return "No eclosion to adulthood"

# --- Scenario 1 ---
# Mothers on 250mg/L cholesterol
# Offspring on 250mg/L cholestanol
s1_maternal_diet = {'cholesterol': 250}
s1_offspring_diet = {'cholestanol': 250, 'cholesterol': 0}
outcome1 = predict_drosophila_outcome(s1_maternal_diet, s1_offspring_diet)

print("--- Scenario 1 Analysis ---")
print(f"Maternal diet: {s1_maternal_diet['cholesterol']}mg/L cholesterol")
print(f"Offspring diet: {s1_offspring_diet['cholestanol']}mg/L cholestanol and {s1_offspring_diet['cholesterol']}mg/L cholesterol")
print(f"Predicted outcome: {outcome1}\n")

# --- Scenario 2 ---
# Mothers on 250mg/L cholestanol
# Offspring on 250mg/L cholestanol and 2mg/L cholesterol
s2_maternal_diet = {'cholestanol': 250, 'cholesterol': 0}
s2_offspring_diet = {'cholestanol': 250, 'cholesterol': 2}
outcome2 = predict_drosophila_outcome(s2_maternal_diet, s2_offspring_diet)

print("--- Scenario 2 Analysis ---")
print(f"Maternal diet: {s2_maternal_diet['cholestanol']}mg/L cholestanol (implies no maternal cholesterol provisioning)")
print(f"Offspring diet: {s2_offspring_diet['cholestanol']}mg/L cholestanol and {s2_offspring_diet['cholesterol']}mg/L cholesterol")
print(f"Predicted outcome: {outcome2}\n")

print(f"Final Answer Pair: ({outcome1}, {outcome2})")