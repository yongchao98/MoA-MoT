import math

def demonstrate_probability_calculations():
    """
    This function demonstrates the flawed calculation made by Agent A
    and contrasts it with a correct probabilistic approach.
    """
    num_quarters = 20
    base_prob_increase = 0.01

    # --- Agent A's Flawed Calculation ---
    print("--- Agent A's Flawed Calculation (Incorrect Linear Addition) ---")
    print("Agent A incorrectly assumed the total probability is the sum of the quarterly probabilities.")
    
    # Building the equation string with each number
    equation_parts_A = [f"{base_prob_increase:.2f}" for _ in range(num_quarters)]
    equation_str_A = " + ".join(equation_parts_A)
    
    # The flawed calculation
    flawed_total_prob = base_prob_increase * num_quarters
    
    print(f"Flawed P(Failure) = {equation_str_A} = {flawed_total_prob:.2f}")
    print("\nThis reasoning is flawed. Probabilities for non-mutually-exclusive events cannot be summed directly.\n")


    # --- A More Correct Calculation (Using A's Own Model) ---
    print("--- A Correct Probabilistic Calculation (Using A's Model) ---")
    print("The correct way is to calculate the total probability of success and subtract from 1.")
    print("We'll use A's model where the conditional risk of failure in quarter 'i' is i * 0.01.")

    prob_total_survival = 1.0
    quarterly_survival_prob_strs = []
    
    # Loop through each quarter to calculate cumulative survival
    for i in range(1, num_quarters + 1):
        # A's model: risk in quarter i is i * 0.01
        prob_failure_in_quarter = i * base_prob_increase
        # The probability of survival in that specific quarter
        prob_survival_in_quarter = 1 - prob_failure_in_quarter
        
        # Add the term to the equation string, handling potential negative probabilities in this extreme model
        if prob_survival_in_quarter >= 0:
            quarterly_survival_prob_strs.append(f"{prob_survival_in_quarter:.2f}")
        else:
            # This case won't be hit with num_quarters=20, but it's good practice
            quarterly_survival_prob_strs.append(f"(1 - {prob_failure_in_quarter:.2f})")
        
        # Update the total survival probability
        prob_total_survival *= prob_survival_in_quarter

    # Build the equation string for total survival
    equation_str_survival = " * ".join(quarterly_survival_prob_strs)
    
    print(f"P(Total Success) = {equation_str_survival}")
    print(f"Calculated P(Total Success) = {prob_total_survival:.4f}\n")

    correct_total_prob_failure = 1 - prob_total_survival
    print("P(Total Failure) = 1 - P(Total Success)")
    print(f"Correct P(Failure) = 1 - {prob_total_survival:.4f} = {correct_total_prob_failure:.4f}")
    print("\nThis shows that even with A's pessimistic model, the actual risk, when calculated correctly, is much higher than his flawed 0.2 estimate.")
    

    # --- Calculation for the Actual Reported Risk ---
    print("\n--- Analysis of Actual Risk (≤ 0.05) ---")
    print("To achieve a total failure risk of 5% over 20 quarters with a constant quarterly risk 'p':")
    
    actual_max_risk = 0.05
    # Solve for p in: 1 - (1 - p)^20 = 0.05
    constant_quarterly_risk = 1 - (1 - actual_max_risk)**(1/num_quarters)
    
    print(f"The required constant quarterly risk 'p' would be 1 - (1 - {actual_max_risk})**(1/{num_quarters})")
    print(f"p ≈ {constant_quarterly_risk:.4f}, or about {constant_quarterly_risk:.2%}.")
    print("This shows the actual quarterly risk was likely very small and did not increase as dramatically as A had assumed.")

demonstrate_probability_calculations()