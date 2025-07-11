import sympy as sp
from fractions import Fraction

def solve_optimal_combination():
    """
    Finds the optimal linear combination of composite Simpson's 1/3 and 3/8 rules
    on an interval [a,b] with 6 subintervals.
    It calculates the constants (C, n, m) for the error term C*(b-a)^n * f^(m)(xi).
    """
    
    # Set up symbolic variables
    # Let the interval be [-d, d], so b-a = 2d
    d = sp.Symbol('d', positive=True)
    x = sp.Symbol('x')
    
    # Define a generic function f by its Taylor series around 0
    # f_coeffs are the derivatives at 0, i.e., f_k = f^(k)(0)
    f_coeffs = [sp.Symbol(f'f{i}') for i in range(10)]
    f_series = sum(c * x**i / sp.factorial(i) for i, c in enumerate(f_coeffs))

    # True Integral over [-d, d]
    I = sp.integrate(f_series, (x, -d, d))
    
    # For 6 subintervals, the step size h is (2d)/6 = d/3
    h = d / 3
    points = [-3*h, -2*h, -h, 0, h, 2*h, 3*h]

    # Composite Simpson's 1/3 rule (3 applications)
    # Weights: 1 4 2 4 2 4 1
    s13_coeffs = [1, 4, 2, 4, 2, 4, 1]
    f_vals = [f_series.subs(x, p) for p in points]
    S13 = (h/3) * sum(c * v for c, v in zip(s13_coeffs, f_vals))

    # Composite Simpson's 3/8 rule (2 applications)
    # Weights: 1 3 3 2 3 3 1
    s38_coeffs = [1, 3, 3, 2, 3, 3, 1]
    S38 = (3*h/8) * sum(c * v for c, v in zip(s38_coeffs, f_vals))
    
    # Error terms E = I - S
    E13 = sp.simplify(I - S13)
    E38 = sp.simplify(I - S38)
    
    # Extract coefficients of the leading error term (d^5 * f4)
    E13_f4_coeff = E13.coeff(f_coeffs[4]).coeff(d**5)
    E38_f4_coeff = E38.coeff(f_coeffs[4]).coeff(d**5)
    
    # Solve for alpha and beta in S_opt = alpha*S13 + beta*S38
    # alpha + beta = 1
    # alpha * E13_coeff + beta * E38_coeff = 0
    alpha, beta = sp.symbols('alpha beta')
    eq1 = sp.Eq(alpha + beta, 1)
    eq2 = sp.Eq(alpha * E13_f4_coeff + beta * E38_f4_coeff, 0)
    solution = sp.solve([eq1, eq2], [alpha, beta])
    
    alpha_val = solution[alpha]
    beta_val = solution[beta]

    # Calculate the total error of the optimal rule
    E_opt = sp.simplify(alpha_val * E13 + beta_val * E38)

    # Find the leading term of the new error
    m = 0
    for i in range(len(f_coeffs)):
        if E_opt.coeff(f_coeffs[i]) != 0:
            m = i
            break
            
    n = 0
    if m > 0:
        error_term = E_opt.coeff(f_coeffs[m])
        # Find the power of d
        n_d = sp.degree(error_term, d)
        # The error term is of the form K * d^n_d * f^(m)
        # Since d = (b-a)/2, d^n_d = ((b-a)/2)^n_d = (b-a)^n_d / 2^n_d
        # So n = n_d
        n = n_d
        
        # Get the constant K
        K_d = error_term.coeff(d**n)
        
        # The question's error form is C * (b-a)^n * f^(m)
        # Our error is K_d * d^n * f^(m) = K_d * ((b-a)/2)^n * f^(m)
        # So C*(b-a)^n = K_d * (b-a)^n / 2^n  => C = K_d / 2^n
        C_val = K_d / (2**n)

    # The standard error definition I-S gives a negative C.
    # The question specifies C>0, implying error is S-I.
    C_final = abs(C_val)

    print(f"The optimal linear combination is S = ({alpha_val})*S_1/3 + ({beta_val})*S_3/8.")
    print(f"The degree of precision is {m-1}.")
    print(f"The error term is of the form C * (b-a)^n * f^(m)(xi).")
    print(f"Calculated values are:")
    print(f"m (order of derivative) = {m}")
    print(f"n (power of b-a) = {n}")
    # Use Fraction for exact representation
    C_frac = Fraction(C_final).limit_denominator()
    print(f"C = {C_frac.numerator}/{C_frac.denominator}")

    # For submission format
    # print(f"({C_frac.numerator}/{C_frac.denominator}, {n}, {m})")

solve_optimal_combination()
<<< (1/39191040, 7, 6) >>>