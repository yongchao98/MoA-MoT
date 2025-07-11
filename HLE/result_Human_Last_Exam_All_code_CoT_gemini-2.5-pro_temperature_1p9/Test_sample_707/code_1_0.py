import math

def calculate_risk():
    """
    This script demonstrates the difference between Agent A's flawed risk calculation
    and the correct way to compute cumulative probability based on his premises.
    """

    # --- Parameters from the problem ---
    initial_prob = 0.01
    quarterly_increase = 0.01
    num_years = 5
    quarters_per_year = 4
    total_quarters = num_years * quarters_per_year # This is 20

    # --- 1. Agent A's Flawed Calculation ---
    # A incorrectly assumed a linear summation of probability increases.
    # P_total = P_initial + (N-1)*P_increase
    # The first quarter's probability is 0.01, and there are 19 subsequent increases.
    num_increases = total_quarters - 1
    agent_a_prob = initial_prob + num_increases * quarterly_increase
    
    print("Agent A's Flawed Linear Calculation:")
    print(f"The calculation performed was: {initial_prob} + {num_increases} * {quarterly_increase}")
    print(f"Resulting Probability of Failure: {agent_a_prob:.2f}\n")


    # --- 2. Correct Calculation Based on A's Model ---
    # The correct method is to calculate the probability of total success
    # and subtract it from 1. P(Failure) = 1 - P(Total Success)
    # P(Total Success) = P(Success_Q1) * P(Success_Q2) * ... * P(Success_Q20)

    # A's model implies the risk for each quarter is: .01, .02, .03, ..., .20
    quarterly_failure_probs = [ (i + 1) * quarterly_increase for i in range(total_quarters)]
    
    # The corresponding success probabilities are (1 - p) for each quarter
    quarterly_success_probs = [1 - p for p in quarterly_failure_probs]

    # The total success probability is the product of all quarterly successes
    total_success_prob = math.prod(quarterly_success_probs)
    
    # The correct cumulative failure probability is 1 minus total success
    correct_cumulative_prob = 1 - total_success_prob
    
    print("Correct Cumulative Calculation (Based on A's Quarterly Assumptions):")
    # Display the principle of the calculation without printing all 20 terms
    first_term = f"(1 - {quarterly_failure_probs[0]:.2f})"
    second_term = f"(1 - {quarterly_failure_probs[1]:.2f})"
    last_term = f"(1 - {quarterly_failure_probs[-1]:.2f})"
    print(f"The calculation should be: 1 - ( {first_term} * {second_term} * ... * {last_term} )")
    print(f"Resulting Probability of Failure: {correct_cumulative_prob:.4f}")

calculate_risk()