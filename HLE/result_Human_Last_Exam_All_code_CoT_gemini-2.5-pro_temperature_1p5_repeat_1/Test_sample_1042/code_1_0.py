import collections

def identify_E_Ya_given_AL(a, a_prime, l):
    """
    Calculates the identifiable value of E(Y^a | A=a', L=l).

    Args:
        a (int): The counterfactual treatment value (0 or 1).
        a_prime (int): The observed treatment value to condition on (0 or 1).
        l (int): The value of the confounder L (0 or 1).
    """

    # --- Pre-computed quantities from a hypothetical dataset ---

    # P(A=a' | L=l)
    # Format: P_A_given_L[l][a']
    P_A_given_L = {
        0: {0: 0.6, 1: 0.4},
        1: {0: 0.3, 1: 0.7}
    }

    # E(Y | A=a', L=l)
    # Format: E_Y_given_AL[l][a']
    E_Y_given_AL = {
        0: {0: 10.0, 1: 15.0},
        1: {0: 12.0, 1: 20.0}
    }

    # E(Y^a | L=l), assumed to be identifiable as per the problem statement
    # Format: E_Ya_given_L[l][a]
    E_Ya_given_L = {
        0: {0: 11.0, 1: 14.0},
        1: {0: 13.0, 1: 18.0}
    }

    # --- Identification Logic ---

    # Case 1: a == a_prime (using Consistency)
    # E(Y^a | A=a, L=l) = E(Y | A=a, L=l)
    if a == a_prime:
        result = E_Y_given_AL[l][a]
        print(f"For a={a}, a'={a_prime}, l={l}:")
        print(f"E(Y^a={a} | A={a_prime}, L={l}) = E(Y | A={a_prime}, L={l}) = {result:.2f}")
        print("-" * 30)
        return

    # Case 2: a != a_prime (using Law of Total Expectation)
    # E(Y^a | A=a',L) = [E(Y^a|L) - P(A=a|L)E(Y^a|A=a,L)] / P(A=a'|L)
    # where E(Y^a|A=a,L) = E(Y|A=a,L) by consistency.
    
    # Check for positivity/overlap
    if P_A_given_L[l][a_prime] == 0:
        print(f"Cannot identify for a={a}, a'={a_prime}, l={l} due to P(A={a_prime}|L={l})=0 (positivity violation).")
        print("-" * 30)
        return
        
    e_ya_given_l = E_Ya_given_L[l][a]
    p_a_given_l = P_A_given_L[l][a]
    e_y_given_al_consistent = E_Y_given_AL[l][a]
    p_a_prime_given_l = P_A_given_L[l][a_prime]

    numerator = e_ya_given_l - p_a_given_l * e_y_given_al_consistent
    denominator = p_a_prime_given_l
    result = numerator / denominator
    
    print(f"For a={a}, a'={a_prime}, l={l}:")
    print(f"E(Y^a={a} | A={a_prime}, L={l}) = [E(Y^a={a}|L={l}) - P(A={a}|L={l}) * E(Y|A={a},L={l})] / P(A={a_prime}|L={l})")
    print(f"                  = ({e_ya_given_l} - {p_a_given_l} * {e_y_given_al_consistent}) / {p_a_prime_given_l}")
    print(f"                  = ({e_ya_given_l} - {p_a_given_l * e_y_given_al_consistent:.2f}) / {p_a_prime_given_l}")
    print(f"                  = {numerator:.2f} / {denominator}")
    print(f"                  = {result:.2f}")
    print("-" * 30)


if __name__ == '__main__':
    print("Identifying E(Y^a | A, L) for all combinations of binary a, A, and L.\n")
    # Iterate through all combinations for l, a, and a_prime
    for l_val in [0, 1]:
        for a_val in [0, 1]:
            for a_prime_val in [0, 1]:
                identify_E_Ya_given_AL(a=a_val, a_prime=a_prime_val, l=l_val)