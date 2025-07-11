import sympy

def find_stabilizing_controller():
    """
    This function outlines the structure of all proper stabilizing controllers H_2(s)
    for the system H_1(s) = s / (s^2 - 1) using the Youla-Kucera parametrization.
    The final answer is presented as a transfer function H_2(s) = N_c(s) / D_c(s)
    parametrized by a stable, proper rational function K(s).
    """
    # Define the symbolic variable 's' and the parameter function 'K(s)'
    s = sympy.Symbol('s')
    K = sympy.Function('K')(s)

    # Based on the Youla-Kucera parametrization method detailed in the plan,
    # the controller H_2(s) is constructed. The steps involve coprime factorization
    # and solving the Bezout identity, leading to the general form:
    # H_2(s) = (X_0 + D*K) / (Y_0 - N*K).
    #
    # With:
    # N(s) = s/(s+1)^2
    # D(s) = (s-1)/(s+1)
    # X_0(s) = 4
    # Y_0(s) = (s-1)/(s+1)
    #
    # Substituting these into the formula and simplifying by multiplying the numerator
    # and denominator by (s+1)^2 gives the final polynomial form.

    # Numerator of H_2(s) after simplification
    num_term_s2 = 4
    num_term_s = 8
    num_term_const = 4
    k_num_term_s2 = 1
    k_num_term_const = -1
    
    # Denominator of H_2(s) after simplification
    den_term_s2 = 1
    den_term_const = -1
    k_den_term_s = -1

    print("The set of all proper stabilizing controllers H_2(s) is parametrized by a stable, proper rational function K(s).")
    print("The controller has the form H_2(s) = N_c(s) / D_c(s), where:\n")

    # Print the numerator with explicit coefficients
    print(f"Numerator: N_c(s) = ({num_term_s2})*s**2 + ({num_term_s})*s + ({num_term_const}) + (({k_num_term_s2})*s**2 + ({k_num_term_const}))*K(s)")

    # Print the denominator with explicit coefficients
    print(f"\nDenominator: D_c(s) = ({den_term_s2})*s**2 + ({den_term_const}) + (({k_den_term_s})*s)*K(s)")


find_stabilizing_controller()
<<<H_2(s) = \frac{4s^2+8s+4 + (s^2-1)K(s)}{s^2-1 - sK(s)}>>>