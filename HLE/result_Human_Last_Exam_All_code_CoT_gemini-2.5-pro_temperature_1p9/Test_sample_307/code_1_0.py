def solve_drosophila_puzzle():
    """
    Analyzes two Drosophila development scenarios based on provided dietary rules.
    """

    # --- Define Scenarios ---
    # Scenario 1:
    maternal_diet_1_cholesterol = 250
    offspring_diet_1_cholestanol = 250

    # Scenario 2:
    offspring_diet_2_cholesterol = 2
    offspring_diet_2_cholestanol = 250
    # The mother's diet is 250mg/L cholestanol.

    # --- Logic ---
    # For Scenario 1:
    # Mothers on 250mg/L cholesterol provide a store in the egg.
    # Offspring diet of 250mg/L cholestanol is unusable.
    # Once the maternal store is depleted, the larvae cannot produce molting hormone and will die.
    outcome_1 = "Death"

    # For Scenario 2:
    # Mothers on cholestanol provide no cholesterol stores.
    # Offspring diet has 2mg/L cholesterol.
    # The problem explicitly states that larvae "cannot survive to adulthood" on 2mg/L cholesterol.
    outcome_2 = "No eclosion to adulthood"


    # --- Print Results ---
    print("Analysis of Drosophila Development")
    print("="*40)

    print("Question 1:")
    print(f"Mothers are reared on {maternal_diet_1_cholesterol}mg/L cholesterol.")
    print(f"Offspring are fed on {offspring_diet_1_cholestanol}mg/L cholestanol.")
    print(f"Predicted Outcome: {outcome_1}")
    print("-"*40)

    print("Question 2:")
    print("Mothers are reared on 250mg/L cholestanol.")
    print(f"Offspring are fed on {offspring_diet_2_cholestanol}mg/L cholestanol and {offspring_diet_2_cholesterol}mg/L cholesterol.")
    print(f"Predicted Outcome: {outcome_2}")
    print("="*40)

solve_drosophila_puzzle()