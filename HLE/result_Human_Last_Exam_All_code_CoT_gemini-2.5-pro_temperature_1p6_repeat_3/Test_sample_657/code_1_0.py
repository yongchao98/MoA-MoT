import sympy

def find_stabilizing_controllers():
    """
    This function determines the set of all proper stabilizing controllers H_2(s)
    for the plant H_1(s) = s / (s^2 - 1) using Youla-Kucera parametrization.
    """
    # Define the symbolic variable 's' and the parameter K(s)
    s = sympy.Symbol('s')
    K = sympy.Function('K')(s)

    # Step 1: Define the plant's numerator and denominator polynomials.
    n_p = s
    d_p = s**2 - 1

    # Step 2: Perform a coprime factorization of H_1(s) into N(s)/D(s)
    # where N and D are stable and proper.
    # We choose a stable polynomial a(s) = (s+2)^2 to cancel the instabilities.
    a = (s + 2)**2
    N = n_p / a
    D = d_p / a

    # Step 3: Solve the Bezout identity X*N + Y*D = 1 for stable, proper X and Y.
    # This is equivalent to s*X + (s^2-1)*Y = (s+2)^2.
    # From the analysis (see plan), we derive a particular solution:
    X = (14*s + 13) / (s + 2)
    Y = (s - 8) / (s + 2)

    # Step 4: Formulate the family of stabilizing controllers H_2(s).
    # H_2(s) = (X + D*K) / (Y - N*K)
    # To simplify, we combine terms over a common denominator (s+2)^2.
    
    num_H2 = ((X * a) + (D * a * K)).simplify()
    den_H2 = ((Y * a) - (N * a * K)).simplify()

    # Expand the polynomials for a clear final expression.
    num_H2_expanded = sympy.expand(num_H2)
    den_H2_expanded = sympy.expand(den_H2)
    
    # Extract the polynomial parts for printing.
    num_poly_part = sympy.poly(num_H2_expanded.subs(K, 0), s)
    num_k_part = sympy.poly(sympy.collect(num_H2_expanded, K).coeff(K), s)
    
    den_poly_part = sympy.poly(den_H2_expanded.subs(K, 0), s)
    den_k_part = sympy.poly(sympy.collect(den_H2_expanded, K).coeff(K), s)
    
    # Print the explanation and the final result.
    print("The set of all proper stabilizing controllers H_2(s) is given by the parametrization:")
    print("H_2(s) = Numerator(s) / Denominator(s)\n")
    print("where K(s) must be an arbitrary stable and proper rational function.")
    print("The final parametrized transfer function is:\n")

    # The following print statement explicitly shows each number in the final equation.
    print(f"Numerator(s) = ({num_poly_part.expr.coeff(s, 2)}*s^2 + {num_poly_part.expr.coeff(s, 1)}*s + {num_poly_part.expr.coeff(s, 0)}) + ({num_k_part.expr.coeff(s, 2)}*s^2 + {num_k_part.expr.coeff(s, 1)}*s + {num_k_part.expr.coeff(s, 0)})*K(s)")
    print(f"Denominator(s) = ({den_poly_part.expr.coeff(s, 2)}*s^2 + {den_poly_part.expr.coeff(s, 1)}*s + {den_poly_part.expr.coeff(s, 0)}) + ({den_k_part.expr.coeff(s, 1)}*s + {den_k_part.expr.coeff(s, 0)})*K(s)")


find_stabilizing_controllers()