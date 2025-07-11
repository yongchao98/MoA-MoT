def get_drosophila_outcome(offspring_diet, maternal_diet):
    """
    Determines the developmental outcome for Drosophila based on diet.

    Args:
        offspring_diet (dict): A dictionary describing the offspring's diet.
                               Keys: 'cholesterol', 'cholestanol'. Values in mg/L.
        maternal_diet (dict): A dictionary describing the mother's diet.
                              Keys: 'cholesterol', 'cholestanol'. Values in mg/L.

    Returns:
        str: The predicted outcome.
    """
    # Rule 1: Cholestanol is not a usable sterol source.
    # The effective cholesterol is only from the cholesterol component of the diet.
    effective_cholesterol = offspring_diet.get('cholesterol', 0)

    # Rule 2: Based on the text for different cholesterol levels.
    if effective_cholesterol >= 250:
        return "Normal development"
    elif effective_cholesterol > 0 and effective_cholesterol <= 2:
        # Text: "larvae cannot survive to adulthood on diets of 2mg/L cholesterol"
        return "No eclosion to adulthood"
    else:  # effective_cholesterol is 0
        # Maternal contribution can allow initial survival but is not sufficient
        # for a full life cycle on a zero-sterol diet.
        # Text: larvae cannot survive on 0mg/L cholesterol.
        # Text: Adult survival is zero on 250mg/L cholestanol diets.
        # This points to a lethal outcome.
        return "Death"


def solve_scenarios():
    """
    Solves the two scenarios presented in the problem.
    """
    # Scenario 1:
    # Mothers on 250mg/L cholesterol.
    # Offspring on 250mg/L cholestanol.
    maternal_diet_1 = {'cholesterol': 250, 'cholestanol': 0}
    offspring_diet_1 = {'cholesterol': 0, 'cholestanol': 250}
    outcome_1 = get_drosophila_outcome(offspring_diet_1, maternal_diet_1)

    print("--- Scenario 1 ---")
    print(f"Maternal Diet: {maternal_diet_1['cholesterol']}mg/L cholesterol")
    print(f"Offspring Diet: {offspring_diet_1['cholesterol']}mg/L cholesterol, {offspring_diet_1['cholestanol']}mg/L cholestanol")
    print(f"Conclusion 1: The diet provides effectively 0mg/L of usable sterol. The larvae will eventually die.")
    print(f"Predicted Outcome 1: {outcome_1}")
    print("-" * 20)

    # Scenario 2:
    # Mothers on 250mg/L cholestanol.
    # Offspring on 250mg/L cholestanol and 2mg/L cholesterol.
    maternal_diet_2 = {'cholesterol': 0, 'cholestanol': 250}
    offspring_diet_2 = {'cholesterol': 2, 'cholestanol': 250}
    outcome_2 = get_drosophila_outcome(offspring_diet_2, maternal_diet_2)

    print("--- Scenario 2 ---")
    print(f"Maternal Diet: {maternal_diet_2['cholestanol']}mg/L cholestanol")
    print(f"Offspring Diet: {offspring_diet_2['cholesterol']}mg/L cholesterol, {offspring_diet_2['cholestanol']}mg/L cholestanol")
    print(f"Conclusion 2: The diet provides only 2mg/L of cholesterol, which the text states is insufficient for survival to adulthood.")
    print(f"Predicted Outcome 2: {outcome_2}")
    print("-" * 20)
    
    print(f"\nFinal Answer Pair: ({outcome_1}, {outcome_2})")


if __name__ == "__main__":
    solve_scenarios()
<<<C>>>