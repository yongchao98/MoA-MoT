import sympy as sp
from fractions import Fraction

def solve_optimal_quadrature():
    """
    Finds the constants C, n, m for the error term of an optimal linear
    combination of Simpson's 1/3 and Simpson's 3/8 rules.
    """
    # Define symbols for Taylor expansion
    h, c = sp.symbols('h c')
    # Use f_i to represent the i-th derivative of f at c, i.e., f^(i)(c)
    f_derivs = [sp.Symbol(f"f_{i}") for i in range(11)]

    # Taylor series for f(x) around c
    def f(x):
        return sum(f_derivs[i] * (x - c)**i / sp.factorial(i) for i in range(11))

    # 1. Expand the true integral I = integral from c-h to c+h of f(x) dx
    I = sp.integrate(f(x), (x, c-h, c+h)).expand()

    # 2. Expand Simpson's 1/3 rule approximation
    # S_1/3 = (b-a)/6 * [f(a) + 4f((a+b)/2) + f(b)]
    # With a=c-h, b=c+h, (b-a)=2h, (a+b)/2=c
    # S_1/3 = 2h/6 * [f(c-h) + 4f(c) + f(c+h)]
    S_1_3 = (sp.sympify(2)*h/6 * (f(c-h) + 4*f(c) + f(c+h))).expand()

    # 3. Expand Simpson's 3/8 rule approximation
    # S_3/8 = (b-a)/8 * [f(a) + 3f(a+(b-a)/3) + 3f(a+2(b-a)/3) + f(b)]
    # With a=c-h, b=c+h, points are c-h, c-h/3, c+h/3, c+h
    S_3_8 = (sp.sympify(2)*h/8 * (f(c-h) + 3*f(c-h/3) + 3*f(c+h/3) + f(c+h))).expand()

    # 4. Calculate the error series E = I - S for each rule
    E_1_3 = (I - S_1_3).simplify()
    E_3_8 = (I - S_3_8).simplify()
    
    # 5. Find the coefficients of the leading error term (h^5*f_4) for both rules
    c5_1_3 = E_1_3.coeff(h**5 * f_derivs[4])
    c5_3_8 = E_3_8.coeff(h**5 * f_derivs[4])
    
    # 6. Solve for alpha in alpha*c5_1_3 + (1-alpha)*c5_3_8 = 0
    # alpha * (c5_1_3 - c5_3_8) = -c5_3_8
    alpha = -c5_3_8 / (c5_1_3 - c5_3_8)
    beta = 1 - alpha

    # 7. Find the coefficient of the next error term (h^7*f_6)
    c7_1_3 = E_1_3.coeff(h**7 * f_derivs[6])
    c7_3_8 = E_3_8.coeff(h**7 * f_derivs[6])
    
    # 8. Calculate the new coefficient for the h^7 term
    new_c7 = (alpha * c7_1_3 + beta * c7_3_8).simplify()
    
    # The new leading error term is new_c7 * h^7 * f_6
    # We need to express this in terms of (b-a). Since h = (b-a)/2:
    # Error = new_c7 * ((b-a)/2)^7 * f_6
    # Error = new_c7 / (2^7) * (b-a)^7 * f_6
    final_coeff = new_c7 / (2**7)

    # The problem asks for the error as C*(b-a)^n * f^(m) with C > 0.
    # Our calculation gives E = I - S. If final_coeff is negative, we
    # can define the error as S - I to make the constant positive.
    if final_coeff < 0:
        C = -final_coeff
    else:
        C = final_coeff
        
    n = 7  # from (b-a)^7
    m = 6  # from f_6, which is f^(6)

    # Convert C to a fraction for exact representation
    C_frac = Fraction(C).limit_denominator()
    
    print(f"The optimal combination is S = ({sp.pretty(alpha)})*S_1/3 + ({sp.pretty(beta)})*S_3/8.")
    print("The error term for this optimal rule is E = C * (b-a)^n * f^(m)(xi).")
    print("The constants are:")
    print(f"C = {C_frac}")
    print(f"n = {n}")
    print(f"m = {m}")

solve_optimal_quadrature()