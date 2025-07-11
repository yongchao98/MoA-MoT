import math

def demonstrate_risk_calculation_flaw():
    """
    This script demonstrates the flaw in Agent A's risk calculation by
    comparing his additive method with the correct multiplicative method for
    compounding probabilities over time.
    """
    
    # Parameters from the problem
    num_quarters = 20
    # Agent A's (overestimated) assessment of risk for each quarter
    assumed_quarterly_risk = 0.01

    # --- 1. Agent A's Flawed Calculation ---
    # Agent A incorrectly assumes that the total probability of failure
    # is the sum of the risks of each quarter.
    
    # A's flawed calculation: risk_q1 + risk_q2 + ... + risk_q20
    # This simplifies to: num_quarters * assumed_quarterly_risk
    agent_a_total_risk = num_quarters * assumed_quarterly_risk
    
    print("Agent A's Flawed Calculation (Additive Method):")
    # To represent the sum 0.01 + 0.01 + ... + 0.01 (20 times), we show it as a multiplication
    print(f"The flawed reasoning treats cumulative probability as a simple sum.")
    print(f"Flawed Equation: {num_quarters} * {assumed_quarterly_risk} = {agent_a_total_risk:.2f}")
    print("-" * 30)

    # --- 2. Correct Probabilistic Calculation ---
    # The correct way is to calculate the probability of survival through all quarters
    # and subtract the result from 1.
    # P(Failure) = 1 - P(Survival)
    # P(Survival) = P(Survive Q1) * P(Survive Q2) * ... * P(Survive Q20)
    # P(Survival) = (1 - risk_q1) * (1 - risk_q2) * ...
    
    prob_survival_per_quarter = 1 - assumed_quarterly_risk
    total_prob_survival = math.pow(prob_survival_per_quarter, num_quarters)
    correct_total_risk = 1 - total_prob_survival

    print("Correct Probabilistic Calculation (Multiplicative Method):")
    print(f"The correct method calculates the total chance of survival and subtracts it from 1.")
    # Printing the full equation with its numbers
    print(f"Correct Equation: 1 - (1 - {assumed_quarterly_risk})^{num_quarters} = {correct_total_risk:.4f}")
    print("-" * 30)
    
    # --- 3. Conclusion ---
    print("Conclusion:")
    print(f"Agent A calculated a risk of {agent_a_total_risk:.2f} (or 20%).")
    print(f"The correct calculation, even using A's own high estimate for quarterly risk, is approximately {correct_total_risk:.3f} (or 18.2%).")
    print("This shows his method was fundamentally flawed. The premise that probabilities can be added this way is incorrect.")


if __name__ == '__main__':
    demonstrate_risk_calculation_flaw()
