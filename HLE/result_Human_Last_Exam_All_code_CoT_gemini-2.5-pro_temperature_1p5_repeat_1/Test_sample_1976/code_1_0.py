import math

def calculate_t_norm_for_odd_n(n):
    """
    Calculates the 1-norm of the correlation matrix T for the state Jn with an odd n.

    This function computes the norm by evaluating the summation formula derived
    from the Bloch decomposition and prints the step-by-step calculation.
    It also prints the result from the simplified analytical formula for verification.
    """
    if n % 2 == 0:
        print(f"Error: The provided n={n} is not odd. This formula is derived for odd n.")
        return

    print(f"Calculating the 1-norm of the correlation matrix T for n = {n}.")
    print("The formula is: ||T||_1 = (1 / (1 + 3^n)) * S")
    print(f"where S = sum_{{n0=0 to {n}}} C(n+1, n0) * |1 + 3^(n0-1)*(-1)^n0| * 3^(n+1-n0)")
    print("-" * 30)

    total_sum_s = 0
    # Loop over n0, the number of identity operators in the Pauli string
    for n0 in range(n + 2):
        # We are summing over alpha != 0, so n0 <= n
        if n0 > n:
            continue
            
        comb_term = math.comb(n + 1, n0)
        
        # We use integer arithmetic to avoid precision issues with 3**-1
        if n0 == 0:
            # Term for n0=0: C(n+1,0) * |1+1/3| * 3^(n+1) = 1 * (4/3) * 3^(n+1) = 4 * 3^n
            abs_val_term_num = 4
            abs_val_term_den = 3
            term_value = 4 * (3**n)
            print(f"Term(n0={n0}): C({n+1},{n0}) * |1 + 1/3| * 3^{n+1} = {comb_term} * ({abs_val_term_num}/{abs_val_term_den}) * {3**(n+1)} = {term_value}")
        else:
            power_of_3 = 3**(n0 - 1)
            sign = (-1)**n0
            if sign == 1: # n0 is even
                abs_val_term = 1 + power_of_3
                op_str = "+"
            else: # n0 is odd
                abs_val_term = abs(1 - power_of_3)
                op_str = "-"
            
            term_value = comb_term * abs_val_term * (3**(n + 1 - n0))
            print(f"Term(n0={n0}): C({n+1},{n0}) * |1 {op_str} 3^{n0-1}| * 3^{n+1-n0} = {comb_term} * {abs_val_term} * {3**(n + 1 - n0)} = {term_value}")
        
        total_sum_s += term_value

    print("-" * 30)
    print(f"The total sum S = {total_sum_s}")
    
    denominator = 1 + 3**n
    final_norm = total_sum_s // denominator
    
    print(f"The normalization factor is 1 / (1 + 3^{n}) = 1 / {denominator}")
    print(f"The 1-norm is ||T||_1 = S / (1 + 3^n) = {total_sum_s} / {denominator} = {final_norm}")
    
    print("\n--- Verification ---")
    analytical_result = 2**(n + 1) - 1
    print(f"The simplified analytical formula for odd n is 2^(n+1) - 1.")
    print(f"For n={n}, this is 2^({n+1}) - 1 = {2**(n+1)} - 1 = {analytical_result}.")
    
    if final_norm == analytical_result:
        print("The calculated result matches the analytical formula.")
    else:
        print("Error: The calculated result does NOT match the analytical formula.")

# You can change the value of n here to any odd integer
# Example usage:
n_value = 5
calculate_t_norm_for_odd_n(n_value)
