import math

def calculate_norm_for_jn(n):
    """
    Calculates the 1-norm of the correlation matrix T for the quantum state J_n
    for a given odd integer n.

    The function computes the norm using the full summation formula and prints
    the intermediate steps, then verifies it with the simplified analytical formula.
    """
    if n % 2 == 0:
        print("This formula is derived for odd n.")
        return

    print(f"Calculating the 1-norm for n = {n}\n")

    # The formula for the 1-norm is:
    # ||T||_1 = (1 / (1 + 3^n)) * sum_{k=1}^{n+1} [C(n+1, k) * 3^k * |1 + (-1)^k * 3^(n-k)|]

    total_sum_of_terms = 0
    
    print("Calculating terms of the sum:")
    # Sum over k from 1 to n+1
    for k in range(1, n + 2):
        # Binomial coefficient C(n+1, k)
        comb = math.comb(n + 1, k)
        
        # The term |1 + (-1)^k * 3^(n-k)|
        abs_term = abs(1 + ((-1)**k) * (3**(n - k)))
        
        # Full term for a given k
        term_value = comb * (3**k) * abs_term
        total_sum_of_terms += term_value
        
        print(f"k = {k}: C({n+1}, {k}) * 3^{k} * |1 + (-1)^{k} * 3^({n}-{k})| = "
              f"{comb} * {3**k} * {abs_term:.4f} = {term_value:.4f}")

    print(f"\nTotal sum of terms = {total_sum_of_terms:.4f}")

    # Normalization factor is 1 + 3^n
    norm_factor = 1 + 3**n
    print(f"Normalization factor = 1 + 3^{n} = {norm_factor}")

    # Final 1-norm
    norm_T1 = total_sum_of_terms / norm_factor
    print(f"\nThe 1-norm is the sum divided by the normalization factor:")
    print(f"||T||_1 = {total_sum_of_terms:.0f} / {norm_factor:.0f} = {norm_T1:.4f}")

    # Verification with the simplified formula: 2^(n+1) - 1
    analytical_result = 2**(n + 1) - 1
    print("\nVerifying with the simplified analytical formula:")
    print(f"Result from formula 2^(n+1) - 1 = 2^({n+1}) - 1 = {analytical_result}")
    
    # Final confirmation
    if math.isclose(norm_T1, analytical_result):
        print("\nThe calculated value matches the analytical result.")
    else:
        print("\nThere is a discrepancy between the calculation and the analytical formula.")

# Let's run the calculation for an odd integer, for example n=5.
calculate_norm_for_jn(5)
