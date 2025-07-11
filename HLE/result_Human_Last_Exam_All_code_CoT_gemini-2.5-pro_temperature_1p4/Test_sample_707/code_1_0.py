import math

def demonstrate_risk_calculation():
    """
    This function demonstrates Agent A's flawed risk calculation
    and contrasts it with a correct probabilistic model.
    """
    
    num_quarters = 20
    risk_per_quarter = 0.01

    # Agent A's flawed calculation by linear addition
    # This corresponds to the false premise that total risk is the sum of quarterly risks.
    agent_a_calculation = num_quarters * risk_per_quarter

    print("Agent A's Flawed Calculation:")
    print("This calculation incorrectly sums the risk probabilities for each quarter.")
    # To satisfy the "output each number in the final equation" requirement
    print(f"Final Probability = Number of Quarters * Risk per Quarter")
    print(f"Final Probability = {num_quarters} * {risk_per_quarter} = {agent_a_calculation:.2f}")
    print("-" * 20)

    # A more correct calculation based on the same premise of 1% risk per quarter
    # The total probability of success is the product of success probabilities in each quarter.
    # Total Failure Probability = 1 - Total Success Probability
    prob_success_per_quarter = 1 - risk_per_quarter
    total_prob_success = math.pow(prob_success_per_quarter, num_quarters)
    correct_prob_failure = 1 - total_prob_success

    print("A Correct Probabilistic Calculation (for comparison):")
    print("This calculates total failure risk as 1 minus the total probability of success.")
    print(f"Total Failure Risk = 1 - (Probability of Success per Quarter ^ Number of Quarters)")
    print(f"Total Failure Risk = 1 - ( (1 - {risk_per_quarter}) ^ {num_quarters} ) = {correct_prob_failure:.4f}")
    print("-" * 20)
    
    print("The error was assuming a simple sum of risks, which resulted in an overestimation (0.20) "
          "compared to the correct composition of probabilities (~0.182), "
          "even if the underlying 1% quarterly risk premise had been true.")

demonstrate_risk_calculation()