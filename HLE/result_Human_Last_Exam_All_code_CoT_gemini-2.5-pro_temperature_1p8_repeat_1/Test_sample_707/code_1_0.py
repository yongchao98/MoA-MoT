def demonstrate_probability_error():
    """
    This function demonstrates the difference between Agent A's flawed probability
    calculation and the correct one.
    """
    
    num_years = 5
    quarters_per_year = 4
    prob_increase_per_quarter = 0.01  # The assumed risk for each quarter
    
    num_quarters = num_years * quarters_per_year
    
    # --- Agent A's flawed calculation (Linear Addition) ---
    # Agent A incorrectly added the probabilities for each quarter together.
    flawed_prob_failure = num_quarters * prob_increase_per_quarter
    
    print("Agent A's Flawed Calculation (Linear Addition):")
    # Using an f-string to show the numbers in the final equation
    print(f"{num_quarters} * {prob_increase_per_quarter} = {flawed_prob_failure:.2f} (or {flawed_prob_failure:.0%})")
    print("-" * 30)

    # --- The correct calculation ---
    # The correct way is to calculate the probability of total success and subtract from 1.
    prob_success_per_quarter = 1 - prob_increase_per_quarter
    total_prob_success = prob_success_per_quarter ** num_quarters
    correct_prob_failure = 1 - total_prob_success

    print("Correct Calculation (Multiplicative Principle):")
    # Using f-strings to show the numbers in the final equation
    print(f"Probability of success per quarter = 1 - {prob_increase_per_quarter} = {prob_success_per_quarter}")
    print(f"Total probability of success = {prob_success_per_quarter}^{num_quarters} = {total_prob_success:.4f}")
    print(f"Correct probability of failure = 1 - {total_prob_success:.4f} = {correct_prob_failure:.4f} (or {correct_prob_failure:.2%})")
    print("\nConclusion: The false premise was assuming probabilities could be linearly added.")

demonstrate_probability_error()