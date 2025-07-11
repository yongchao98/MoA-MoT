def find_root_cause():
    """
    This function analyzes the causal chain of the 2008 RMBS crisis
    to determine the most likely root cause factor.
    """
    print("Analyzing the causal chain for the drop in non-agency RMBS value (2004-2008):")

    # Step 1: The immediate cause of value loss is the failure of underlying assets.
    symptom = "E. Default rates"
    print(f"1. The direct cause of value destruction was soaring '{symptom}'. This is a symptom, not the root cause.")

    # Step 2: High defaults were caused by the poor quality of the loans themselves.
    loan_characteristics = ["C. Average FICO scores (low)", "B. The percent of floating rate debt (high)"]
    print(f"2. High defaults were caused by poor loan characteristics, such as {loan_characteristics[0]} and {loan_characteristics[1]}.")

    # Step 3: The poor quality loans were created by the institutions that issued them.
    root_cause = "F. The quality of the loan issuer and RMBS originator"
    print(f"3. These poor-quality loans were systematically created and sold by the institutions responsible for underwriting.")
    print(f"   Therefore, the root cause was '{root_cause}'.")

    # Final "equation" representing the causal flow.
    # The originator's quality determines the loan characteristics, which in turn drive defaults.
    # The final print statement is formatted to follow the prompt's instructions.
    print("\nFinal Conclusion:")
    final_equation = f"The root cause factor is '{root_cause[0]}'"
    print(final_equation)


find_root_cause()