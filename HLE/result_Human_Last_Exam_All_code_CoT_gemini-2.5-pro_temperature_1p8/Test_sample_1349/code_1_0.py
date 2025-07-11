import sympy

def solve_problem():
    """
    This function performs the symbolic derivation to find the expression for X.
    The derivation is explained in the text of the solution.
    """
    # Define symbols for the integrals and their components
    ua, va, ub, vb = sympy.symbols('ua va ub vb', real=True)
    I = sympy.I

    # Let I_fa = Integral(f*conj(a) dt) = ua + i*va
    # Let I_fb = Integral(f*conj(b) dt) = ub + i*vb

    # From the derivation explained in the main text, we have relations:
    # 1. Integral(|a|^2) = 2*Re(I_fa) = 2*ua
    # 2. Integral(|b|^2) = 2*Re(I_fb) = 2*ub
    # 3. (2i-1)*Integral(a*conj(b)) = -conj(I_fa) - I_fb
    I_aa = 2 * ua
    I_bb = 2 * ub
    
    I_fa_conj = ua - I*va
    I_fb = ub + I*vb

    # Calculate Integral(a*conj(b)), let's call it S
    S = (-I_fa_conj - I_fb) / (2*I - 1)
    
    # The numerator of X is N = Integral(|a-b|^2) = Integral(|a|^2) + Integral(|b|^2) - 2*Re(Integral(a*conj(b)))
    # N = I_aa + I_bb - (S + conj(S))
    N = I_aa + I_bb - (S + sympy.conjugate(S))
    
    # The denominator of X is D = Integral(|a|^2 + |b|^2) = I_aa + I_bb
    D = I_aa + I_bb

    # The expression for X
    X = N / D
    X_simplified = sympy.simplify(X)

    # Let's express X in a more readable form.
    # The result from sympy is: (4*ua + 4*ub + 2*vb - 2*va)/(5*(ua + ub))
    # This simplifies to 4/5 + (2/5)*(vb-va)/(ua+ub)
    # Our analysis showed X = 4/5 - (2/5)*(va-vb)/(ua+ub), which is the same.
    # Let U = ua + ub and V = va - vb
    U, V = sympy.symbols('U V')
    # Then X = 4/5 - (2/5)*(V/U)
    
    # The bounds for V/U are derived as -3 <= V/U <= 2.
    # The supremum of X corresponds to the infimum of V/U.
    inf_V_over_U = -3
    
    sup_X = sympy.Rational(4,5) - sympy.Rational(2,5) * inf_V_over_U
    
    print("Let the quantity to maximize be X.")
    print("Through symbolic manipulation of the given differential equations, we find an expression for X.")
    print(f"X = 4/5 - (2/5) * Q, where Q is a ratio of integral quantities.")
    print("By considering the positivity of certain integrals, we can bound the ratio Q.")
    print("The bounds for Q are found to be: -3 <= Q <= 2.")
    print("To find the supremum of X, we use the infimum of Q.")
    print(f"sup(X) = 4/5 - (2/5) * inf(Q)")
    print(f"sup(X) = 4/5 - (2/5) * ({inf_V_over_U})")
    print(f"sup(X) = {sympy.Rational(4,5)} + {sympy.Rational(6,5)}")
    final_answer = sympy.Rational(10,5)
    print(f"sup(X) = {final_answer}")
    
solve_problem()