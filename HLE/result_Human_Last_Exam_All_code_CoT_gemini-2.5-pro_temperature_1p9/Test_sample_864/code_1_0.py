import sympy

def solve_and_print():
    """
    Based on the derivation, this function prints the final expression for the whitening filter W(D).
    We identified the causal, minimum-phase factor of the channel as F(D) = 1 + (6/5)D + (2/3)D^2.
    The whitening filter that ensures a causal equivalent channel is W(D) = 1 / F(1/D).
    """
    
    # Define the symbol D for the D-transform
    D = sympy.symbols('D')
    
    # Coefficients of the derived causal minimum-phase factor F(D)
    f0 = 1
    f1 = 6/5
    f2 = 2/3
    
    # Construct the expression for F(1/D)
    F_invD = f0 + f1 * (D**-1) + f2 * (D**-2)
    
    # The whitening filter is W(D) = 1 / F(1/D)
    W_D = 1 / F_invD
    
    # The problem asks to output each number in the final equation.
    # The final equation is W(D) = 1 / (f0 + f1*D^-1 + f2*D^-2)
    
    print("The derived whitening filter is W(D) = 1 / F(1/D), where F(D) is the causal, minimum-phase spectral factor of the channel.")
    print("Based on the coefficients q1=2 and q2=2/3, we found F(D) = 1 + (6/5)D + (2/3)D^2.")
    print("\nTherefore, the whitening filter W(D) is:")
    
    # To print the numbers in the final equation explicitly:
    f0_s = sympy.pretty(sympy.Rational(f0))
    f1_s = sympy.pretty(sympy.Rational(f1))
    f2_s = sympy.pretty(sympy.Rational(f2))

    print(f"      1")
    print(f"W(D) = -----------------------------")
    print(f"       {f0_s} + {f1_s}*D^-1 + {f2_s}*D^-2")
    
    print("\nIn equation form:")
    sympy.pprint(W_D, use_unicode=False)

solve_and_print()