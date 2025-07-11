import sympy

def solve():
    """
    This function derives and prints the formula for P(n).
    
    The problem is to find P(n) such that the approximation:
    T_refined(n) = A * n^(L/2) * (1 + L/(2*n) + P(n))
    has a relative error of O((L/n)^4) for Q(n) = Product_{k=1 to n} k^(1/k), where L = ln(n).
    
    This is equivalent to making the difference between their logarithms, 
    ln(Q(n)) - ln(T_refined(n)), of order O((L/n)^4).
    
    Step 1: Asymptotic expansion for ln(Q(n)).
    Using the Euler-Maclaurin formula for f(x) = ln(x)/x, we get:
    ln(Q(n)) ~ C + (L^2)/2 + L/(2*n) + (L-1)/(12*n^2) + (11-6*L)/(720*n^4) + ...
    Let's call the part after C + (L^2)/2 as the "Target Series".
    
    Step 2: Taylor expansion for ln(T_refined(n)).
    ln(T_refined(n)) = ln(A) + (L^2)/2 + ln(1 + L/(2*n) + P(n))
    The constant C from Euler-Maclaurin corresponds to ln(A).
    Let X = L/(2*n) + P(n). The expansion is ln(1+X) = X - X^2/2 + X^3/3 - ...
    
    Step 3: Solve for P(n).
    We set ln(1 + L/(2*n) + P(n)) equal to the Target Series and solve for P(n).
    P(n) can be solved by equating coefficients for each power of 1/n.
    Let P(n) = P2/n^2 + P3/n^3 + ...
    
    For the n^-2 term:
    From ln(1+X): P2 - (L^2)/8
    From Target Series: (L-1)/12
    Equating gives: P2 = (L-1)/12 + (L^2)/8 = (2*L-2 + 3*L^2)/24
    
    For the n^-3 term:
    From ln(1+X): P3 - (L*P2_coeffs)/(2) + (L^3)/24 
    (after discovering a mistake in manual derivation, the coefficient is L/2 not L)
    From Target Series: 0
    Equating gives: P3 = (L*(3*L^2+2*L-2))/(2*24) - L^3/24 = (3*L^3+2*L^2-2*L - 2*L^3)/48 = (L^3+2*L^2-2*L)/48

    By including P2 and P3, the log-error starts at the n^-4 term. The leading part of this error is of order L^4/n^4. So the error is O((L/n)^4).
    """

    # Coefficients for the numerator polynomial of the n^-2 term.
    # Polynomial: c22*L^2 + c21*L + c20
    c22 = 3
    c21 = 2
    c20 = -2
    
    # Denominator for the n^-2 term.
    d2 = 24
    
    # Coefficients for the numerator polynomial of the n^-3 term.
    # Polynomial: c33*L^3 + c32*L^2 + c31*L
    c33 = 1
    c32 = 2
    c31 = -2
    
    # Denominator for the n^-3 term.
    d3 = 48
    
    # Using Sympy for pretty printing the formula
    L, n = sympy.symbols('L n')
    
    p2_num = c22*L**2 + c21*L + c20
    p2_term = p2_num / (d2 * n**2)
    
    p3_num = c33*L**3 + c32*L**2 + c31*L
    p3_term = p3_num / (d3 * n**3)
    
    p_n_expr = p2_term + p3_term
    
    # Print the formula for P(n)
    print(f"P(n) = {sympy.pretty(p_n_expr, use_unicode=False)}")
    # The previous line prints a multiline pretty formula. 
    # Let's also print a compact one-line version for clarity.
    
    term2_str = f"({c22}*L^2 + {c21}*L - {abs(c20)})/({d2}*n^2)"
    term3_str = f"({c33}*L^3 + {c32}*L^2 - {abs(c31)}*L)/({d3}*n^3)"
    
    print("\nOne-line formula:")
    print(f"P(n) = {term2_str} + {term3_str}")

solve()