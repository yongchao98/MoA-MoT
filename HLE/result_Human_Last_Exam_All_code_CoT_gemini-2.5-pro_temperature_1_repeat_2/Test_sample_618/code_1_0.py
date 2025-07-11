import sympy

def display_transformations():
    """
    This function calculates and displays the finite transformation for the x variable
    of the Lie group that leaves the given heat equation invariant.
    """
    # Define symbolic variables
    t, x, epsilon = sympy.symbols('t x epsilon')
    c1, c2, c3, k1 = sympy.symbols('c1 c2 c3 k1', real=True, constant=True)

    # The infinitesimals tau and xi are derived from Lie symmetry analysis
    # tau = c1
    # xi = (c2/k1) * exp(k1*t) + c3

    # We find the finite transformation by integrating the flow equations:
    # dt*/d(epsilon) = tau
    # dx*/d(epsilon) = xi

    print("The general representation for the transformation on x, denoted as x*, is derived from its infinitesimal generator xi.")
    print("The form of the transformation depends on the constant c1 from the time transformation t* = t + c1*epsilon.\n")

    # Case 1: c1 is not equal to 0
    print("Case 1: c1 != 0")
    print("t* = t + c1*epsilon")
    # Integrating dx*/d(epsilon) = (c2/k1)*exp(k1*t*) + c3, where t* = t + c1*epsilon
    # x_star = x + c3*epsilon + (c2/(c1*k1**2)) * (exp(k1*(t + c1*epsilon)) - exp(k1*t))
    t_star_case1 = t + c1 * epsilon
    x_star_case1 = x + c3*epsilon + (c2 / (c1 * k1**2)) * (sympy.exp(k1 * t_star_case1) - sympy.exp(k1*t))
    
    print("The transformation on x is:")
    sympy.pprint(sympy.Eq(sympy.Symbol('x*'), x_star_case1), use_unicode=True)
    # Output each "number" (symbolic constant) in the equation
    print("The constants in this equation are: c1, c2, c3, k1.")
    print("-" * 30)

    # Case 2: c1 = 0
    print("Case 2: c1 = 0")
    print("t* = t")
    # Integrating dx*/d(epsilon) = (c2/k1)*exp(k1*t) + c3
    # x_star = x + ((c2/k1)*exp(k1*t) + c3)*epsilon
    x_star_case2 = x + ((c2/k1)*sympy.exp(k1*t) + c3)*epsilon

    print("The transformation on x is:")
    sympy.pprint(sympy.Eq(sympy.Symbol('x*'), x_star_case2), use_unicode=True)
    # Output each "number" (symbolic constant) in the equation
    print("The constants in this equation are: c2, c3, k1.")

if __name__ == '__main__':
    display_transformations()