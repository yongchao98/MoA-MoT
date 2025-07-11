import sympy

def solve_whitening_filter():
    """
    This function calculates the whitening filter W(D) based on the problem description.
    It assumes the autocorrelation Q(D) arises from a 2nd-degree causal FIR filter G(D),
    and resolves an inconsistency in the provided q_k coefficients by prioritizing
    q_1 and q_2 over q_0.
    """
    
    # Define D as a symbolic variable for display purposes
    D = sympy.Symbol('D')
    
    # Given coefficients from the problem formula for k=1 and k=2
    q1 = 2
    q2 = sympy.Rational(2, 3)
    
    # Assume G(D) = c0 + c1*D + c2*D^2 and normalize with c0=1
    c0 = 1
    
    # From q2 = c0*c2
    c2 = q2 / c0
    
    # From q1 = c1*(c0+c2)
    c1 = q1 / (c0 + c2)
    
    # The causal, minimum-phase factor is G(D) = c0 + c1*D + c2*D^2
    # The roots of 1 + (6/5)D + (2/3)D^2 = 0 are outside the unit circle,
    # so G(D) is minimum-phase.
    
    # The whitening filter that makes the resulting channel causal is W(D) = 1/G(1/D)
    # G(1/D) = c0 + c1*D**-1 + c2*D**-2
    
    print("Based on the analysis, the causal factor G(D) of the channel autocorrelation is found to be:")
    print(f"G(D) = {c0} + ({c1})*D + ({c2})*D^2")
    print("\nThe whitening filter W(D) that ensures the resulting channel is causal is W(D) = 1/G(1/D).")
    print("The final expression for the whitening filter is:")
    
    # We print each number in the final equation as requested.
    # Using f-string for clear output.
    c0_val = c0
    c1_val_frac = sympy.fraction(c1)
    c2_val_frac = sympy.fraction(c2)
    
    print(f"\nW(D) = 1 / ( {c0_val} + ({c1_val_frac[0]}/{c1_val_frac[1]})*D^-1 + ({c2_val_frac[0]}/{c2_val_frac[1]})*D^-2 )")

solve_whitening_filter()
