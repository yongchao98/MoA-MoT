def get_drosophila_outcome(cholesterol_mg_L, cholestanol_mg_L, maternal_health):
    """
    Determines the developmental outcome for Drosophila based on diet.

    Args:
        cholesterol_mg_L (int): Concentration of cholesterol in the diet.
        cholestanol_mg_L (int): Concentration of cholestanol in the diet.
        maternal_health (str): Describes the mother's diet ('Healthy' or 'Deficient').

    Returns:
        str: A description of the developmental outcome.
    """
    print(f"Analyzing diet: {cholesterol_mg_L}mg/L Cholesterol, {cholestanol_mg_L}mg/L Cholestanol. Maternal Health: {maternal_health}")

    # The text states cholestanol is not a usable sterol precursor.
    # Therefore, development depends entirely on the cholesterol concentration.

    if cholesterol_mg_L >= 250:
        return "Normal development"
    elif cholesterol_mg_L == 2:
        # Text explicitly states: "larvae cannot survive to adulthood on diets of 2mg/L cholesterol."
        return "No eclosion to adulthood"
    elif cholesterol_mg_L == 0:
        # A diet with only cholestanol is effectively a 0mg/L cholesterol diet.
        # This is lethal and results in developmental arrest and death.
        return "Death"
    else:
        return "Impossible to determine"

# --- Scenario 1 ---
print("--- Evaluating Scenario 1 ---")
# Mothers on 250mg/L cholesterol are healthy.
# Offspring are on a diet of 250mg/L cholestanol (which means 0mg/L cholesterol).
scenario_1_cholesterol = 0
scenario_1_cholestanol = 250
outcome_1 = get_drosophila_outcome(scenario_1_cholesterol, scenario_1_cholestanol, "Healthy")
print(f"Reasoning: Offspring will use up maternal cholesterol stores from the egg, then rely on their diet.")
print(f"The diet contains {scenario_1_cholesterol}mg/L of usable cholesterol. This is insufficient for molting and leads to developmental arrest.")
print(f"Outcome 1: {outcome_1}\n")


# --- Scenario 2 ---
print("--- Evaluating Scenario 2 ---")
# Mothers on 250mg/L cholestanol are deficient.
# Offspring are on a diet of 2mg/L cholesterol and 250mg/L cholestanol.
scenario_2_cholesterol = 2
scenario_2_cholestanol = 250
outcome_2 = get_drosophila_outcome(scenario_2_cholesterol, scenario_2_cholestanol, "Deficient")
print(f"Reasoning: The deciding factor is the concentration of usable cholesterol, which is {scenario_2_cholesterol}mg/L.")
print(f"The problem states this diet level is not sufficient to reach the adult stage.")
print(f"Outcome 2: {outcome_2}\n")

# --- Final Answer ---
print("--- Conclusion ---")
print(f"The combined outcomes are: '{outcome_1}' and '{outcome_2}'.")
print("This corresponds to answer choice C.")

<<<C>>>