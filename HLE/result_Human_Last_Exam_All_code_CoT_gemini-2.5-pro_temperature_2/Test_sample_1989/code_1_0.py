import sympy

def solve_and_print():
    """
    This script calculates and displays the asymptotic behavior of the solution omega.

    The derivation finds that the solution has the form:
    omega(r, theta) ~ r**(power) * exp(exponent)

    This script constructs and prints this form symbolically.
    """
    
    # 1. Define all symbolic variables
    # r is the radial coordinate, theta is the angle. A and B are the given real constants.
    r, theta, A, B = sympy.symbols('r, theta, A, B', real=True)
    r = sympy.Symbol('r', real=True, positive=True)

    # 2. Define the components of the solution's asymptotic form.
    
    # The exponential decay term, which is the same as in the case A=B=0.
    # It is derived from the constant part of the advection field u = e_1.
    # The term is: exp(-r * (1 - cos(theta)))
    exp_term = sympy.exp(-r * (1 - sympy.cos(theta)))
    
    # The algebraic term is a power of r. The total exponent has several parts.
    
    # a) The original exponent from the A=B=0 case, which is -1/2.
    # This corresponds to the asymptotic behavior of the modified Bessel function K_0.
    base_r_exponent = sympy.Rational(-1, 2)
    
    # b) The correction to the exponent due to the non-zero A and B terms.
    # This part is derived from the WKB analysis of the transformed PDE.
    corrector_exponent = A * (1 - sympy.cos(theta)) - B * sympy.sin(theta)
    
    # The total exponent for r is the sum of the base exponent and the correction.
    total_r_exponent = base_r_exponent + corrector_exponent
    
    # 3. Assemble the full asymptotic expression for omega.
    # omega is proportional to r**(total_exponent) * exp_term
    omega_asymptotic = r**total_r_exponent * exp_term
    
    # 4. The corrector is the factor that multiplies the A=B=0 asymptotic form.
    # The A=B=0 form is r**(-1/2) * exp_term.
    # The corrector is therefore r**(corrector_exponent).
    corrector = r**corrector_exponent
    
    # 5. Print the results in a readable format.
    # sympy.pprint renders mathematical expressions nicely in a terminal.
    
    print("The full large-distance asymptotic behavior for omega(r, theta) is proportional to:")
    # We use use_unicode=False for maximum compatibility with terminals.
    sympy.pprint(omega_asymptotic, use_unicode=False)
    
    print("\n" + "="*50)
    print("The corrector term that multiplies the A=B=0 behavior is:")
    sympy.pprint(corrector, use_unicode=False)

if __name__ == '__main__':
    solve_and_print()