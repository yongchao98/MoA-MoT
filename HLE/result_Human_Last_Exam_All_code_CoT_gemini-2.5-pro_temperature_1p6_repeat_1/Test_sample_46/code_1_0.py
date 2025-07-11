def analyze_rmbs_crisis_root_cause():
    """
    This function analyzes the factors contributing to the decline in non-agency
    RMBS value from 2004-2008 to identify the root cause.
    """

    print("Analyzing the causal chain of the non-agency RMBS crisis:")
    print("----------------------------------------------------------")

    # While all factors played a role, we trace them back to the origin.
    # The problem starts with the creation of low-quality loans.
    origin_of_problem = "The quality of the loan issuer and RMBS originator"
    symptom_1 = "Average FICO scores on the loans were low"
    symptom_2 = "The percent of floating rate debt in the pool was high"
    direct_cause_of_loss = "Default rates soared"
    enabling_factor = "The rating of the credit agency on issuance was inflated"

    print(f"The ultimate root cause was: '{origin_of_problem}'.")
    print("This primary failure led to the creation of poor-quality assets.")
    print("\nConceptual Equation of the Causal Flow:")

    # This 'equation' uses numbers to represent the steps in the causal chain.
    print("Step 1: Poor standards by loan issuers and RMBS originators.")
    print("Step 2: Led to pools of loans with low FICO scores and risky terms.")
    print("Step 3: Which directly resulted in unprecedentedly high default rates.")
    print("Step 4: Causing a total collapse in RMBS values.")
    print("\nFinal Conclusion:")
    print(f"The root cause factor that determined the creation of risky assets and their eventual failure was F: '{origin_of_problem}'.")


if __name__ == "__main__":
    analyze_rmbs_crisis_root_cause()
<<<F>>>