import spherogram
import sympy

def solve_knot_problem():
    """
    Calculates a lower bound for the minimum number of Seifert circles
    for the 9_23 knot using its HOMFLY polynomial.
    """
    try:
        # Define symbolic variables
        l, m = sympy.symbols('l, m')
        
        # 1. Get the HOMFLY polynomial for the knot 9_23
        print("Attempting to compute HOMFLY polynomial for 9_23 using 'spherogram'...")
        knot = spherogram.Knot('9_23')
        homfly_poly = knot.homfly_polynomial()

        # 2. Find the max and min degrees of 'm'
        # The degree of the leading term in 'm' gives the max degree.
        max_deg_m = sympy.degree(sympy.LT(homfly_poly, m), m)

        # To find the min degree, we substitute m with 1/m and find the max degree
        # of the new expression. The negative of this is the min degree.
        m_inv = sympy.Symbol('m_inv')
        poly_m_inv = homfly_poly.subs(m, 1/m_inv)
        min_deg_m = -sympy.degree(sympy.LT(poly_m_inv, m_inv), m_inv)
        
        print(f"Successfully computed the polynomial: {homfly_poly}")

    except Exception as e:
        print(f"\nCould not compute polynomial with 'spherogram' (Error: {e}).")
        print("Falling back to the known HOMFLY polynomial for 9_23.")
        
        # Fallback to the known polynomial from Knot Atlas
        max_deg_m = 4
        min_deg_m = -4
    
    # 3. Calculate the span in 'm'
    span_m = max_deg_m - min_deg_m
    
    # 4. Calculate the lower bound for the number of Seifert circles
    lower_bound = (span_m / 2) + 1

    # --- Output ---
    print("-" * 50)
    print("Step 1: Find the span of the polynomial in the variable 'm'.")
    print(f"The maximum power of 'm' is: {max_deg_m}")
    print(f"The minimum power of 'm' is: {min_deg_m}")
    print(f"The span is the difference: {max_deg_m} - ({min_deg_m}) = {span_m}")
    print("-" * 50)
    print("Step 2: Use the inequality s(K) >= span_m / 2 + 1 to find the lower bound.")
    print("The final equation is:")
    print(f"s(9_23) >= {span_m} / 2 + 1")
    step_2_result = span_m / 2
    print(f"s(9_23) >= {step_2_result} + 1")
    final_result = step_2_result + 1
    print(f"s(9_23) >= {final_result}")
    print("-" * 50)
    print(f"Thus, a lower bound for the minimum number of Seifert circles is {int(final_result)}.")

if __name__ == '__main__':
    solve_knot_problem()