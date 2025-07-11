def get_drosophila_outcome(offspring_diet_cholesterol, offspring_diet_cholestanol, has_maternal_cholesterol):
    """
    Determines the developmental outcome of Drosophila based on diet and maternal contribution.

    Args:
        offspring_diet_cholesterol (int): The mg/L of cholesterol in the offspring's diet.
        offspring_diet_cholestanol (int): The mg/L of cholestanol in the offspring's diet.
        has_maternal_cholesterol (bool): True if eggs have maternally-provided cholesterol.

    Returns:
        str: The predicted outcome.
    """
    # Rule 1: Sufficient cholesterol leads to normal development.
    if offspring_diet_cholesterol >= 250:
        return "Normal development"

    # Rule 2: Insufficient cholesterol (2mg/L) prevents survival to adulthood.
    # This happens regardless of maternal stores, which are not enough for the entire larval stage.
    if offspring_diet_cholesterol == 2:
        return "No eclosion to adulthood"

    # Rule 3: A cholestanol-only or zero-cholesterol diet is lethal.
    # If there's no cholesterol in the diet at all, the larva depends on maternal stores.
    # Once those stores run out, it cannot molt and will die.
    if offspring_diet_cholesterol == 0:
        # The presence of cholestanol is irrelevant as it's not a usable sterol.
        return "Death"

    # Fallback for unexpected cases
    return "Impossible to determine"

# --- Analysis of Scenario 1 ---
print("--- Scenario 1 Analysis ---")
# Mothers were reared on 250mg/L cholesterol.
maternal_diet_1_cholesterol = 250
# Offspring are fed a diet of 250mg/L cholestanol.
offspring_diet_1_cholesterol = 0
offspring_diet_1_cholestanol = 250
# Mothers on a normal diet can provision their eggs.
maternal_stores_1 = True
print(f"Maternal Diet: {maternal_diet_1_cholesterol}mg/L cholesterol.")
print(f"Offspring Diet: {offspring_diet_1_cholesterol}mg/L cholesterol and {offspring_diet_1_cholestanol}mg/L cholestanol.")
outcome_1 = get_drosophila_outcome(offspring_diet_1_cholesterol, offspring_diet_1_cholestanol, maternal_stores_1)
print(f"Predicted Outcome 1: {outcome_1}\n")


# --- Analysis of Scenario 2 ---
print("--- Scenario 2 Analysis ---")
# Mothers were reared on 250mg/L cholestanol.
maternal_diet_2_cholestanol = 250
# Offspring are fed a diet of 250mg/L cholestanol and 2mg/L cholesterol.
offspring_diet_2_cholesterol = 2
offspring_diet_2_cholestanol = 250
# Mothers on a cholestanol diet cannot provide cholesterol to eggs.
maternal_stores_2 = False
print(f"Maternal Diet: {maternal_diet_2_cholestanol}mg/L cholestanol.")
print(f"Offspring Diet: {offspring_diet_2_cholesterol}mg/L cholesterol and {offspring_diet_2_cholestanol}mg/L cholestanol.")
outcome_2 = get_drosophila_outcome(offspring_diet_2_cholesterol, offspring_diet_2_cholestanol, maternal_stores_2)
print(f"Predicted Outcome 2: {outcome_2}\n")

# --- Final Conclusion ---
print(f"The combined outcomes are: ({outcome_1}, {outcome_2})")
print("This corresponds to answer choice C.")

<<<C>>>