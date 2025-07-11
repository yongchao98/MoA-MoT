def get_dietary_rules():
    """Returns a dictionary of rules based on the problem description."""
    rules = {
        "cholesterol_250": "Normal development",
        "cholesterol_2": "No eclosion to adulthood",
        "cholesterol_0": "No eclosion to adulthood",
        "cholestanol_250": "Death" # "Adult survival is zero" implies a lethal outcome.
    }
    return rules

def has_maternal_contribution(maternal_cholesterol_level):
    """Checks if the mother's diet allows for cholesterol contribution to eggs."""
    # A diet of 250mg/L cholesterol is normal and allows for contribution.
    return maternal_cholesterol_level >= 250

def predict_drosophila_outcome(maternal_cholesterol, offspring_cholesterol, offspring_cholestanol):
    """Predicts the developmental outcome based on diets."""
    rules = get_dietary_rules()

    # The long-term outcome is determined by the offspring's own diet,
    # as any maternal contribution is temporary.

    # Check if the diet is explicitly lethal due to cholestanol
    if offspring_cholesterol == 0 and offspring_cholestanol > 0:
        return rules["cholestanol_250"]

    # Check the outcome based on the cholesterol level in the diet
    if offspring_cholesterol == 2:
        return rules["cholesterol_2"]
    elif offspring_cholesterol >= 250:
        return rules["cholesterol_250"]
    else:
        # Covers the 0mg/L case or other undefined insufficient amounts
        return rules["cholesterol_0"]

# --- Scenario 1 ---
mother_diet_1_cholesterol = 250
offspring_diet_1_cholesterol = 0
offspring_diet_1_cholestanol = 250
outcome_1 = predict_drosophila_outcome(mother_diet_1_cholesterol, offspring_diet_1_cholesterol, offspring_diet_1_cholestanol)

print("--- Scenario 1 ---")
print(f"Mothers reared on {mother_diet_1_cholesterol}mg/L cholesterol.")
print(f"Offspring reared on {offspring_diet_1_cholesterol}mg/L cholesterol and {offspring_diet_1_cholestanol}mg/L cholestanol.")
print(f"Predicted Outcome: {outcome_1}")
print("")

# --- Scenario 2 ---
# Mothers on cholestanol have no cholesterol to contribute.
mother_diet_2_cholesterol = 0
offspring_diet_2_cholesterol = 2
offspring_diet_2_cholestanol = 250
outcome_2 = predict_drosophila_outcome(mother_diet_2_cholesterol, offspring_diet_2_cholesterol, offspring_diet_2_cholestanol)

print("--- Scenario 2 ---")
print(f"Mothers reared on {offspring_diet_1_cholestanol}mg/L cholestanol.")
print(f"Offspring reared on {offspring_diet_2_cholesterol}mg/L cholesterol and {offspring_diet_2_cholestanol}mg/L cholestanol.")
print(f"Predicted Outcome: {outcome_2}")
print("")

print(f"The combined answer is: ({outcome_1}, {outcome_2})")