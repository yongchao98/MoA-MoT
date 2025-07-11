import sympy

def find_transformation():
    """
    This function uses the sympy library to define and display the
    finite transformation for the spatial coordinate x, based on the
    Lie group analysis of the given PDE.
    """
    # Define symbols for the variables, constants, and the group parameter
    t, x = sympy.symbols('t x')
    epsilon, c1, c3, c5, k1 = sympy.symbols('epsilon c1 c3 c5 k1')

    # Define the transformed coordinate x* as a function for pretty printing
    x_star = sympy.Symbol('x*')
    
    # The time transformation is t* = t + c1*epsilon
    t_star_expr = t + c1 * epsilon

    # --- Case 1: c1 is not equal to 0 ---
    # The transformation is found by integrating dx*/d(epsilon) = xi(t*)
    # xi(t*) = c3 - (2*c5/k1)*exp(k1*t*)
    # x*(epsilon) = x + integral from 0 to epsilon of xi(t + c1*s) ds
    x_star_expr_c1_nonzero = x + c3 * epsilon - (2 * c5 / (k1**2 * c1)) * (sympy.exp(k1 * t_star_expr) - sympy.exp(k1 * t))

    # --- Case 2: c1 = 0 ---
    # In this case, t* = t, so the integration is simpler
    # x*(epsilon) = x + integral from 0 to epsilon of xi(t) ds
    x_star_expr_c1_zero = x + epsilon * (c3 - (2 * c5 / k1) * sympy.exp(k1 * t))

    # Print the results in a readable format
    print("The general representation for the transformation on the spatial coordinate x, denoted by x*, is derived from the infinitesimal xi.")
    print(f"The infinitesimal is xi = {c3 - (2 * c5 / k1) * sympy.exp(k1 * t)}.")
    print("The finite transformation x* depends on the group parameter epsilon and has two forms, depending on the constant c1 from the time transformation t* = t + c1*epsilon.")
    
    print("\n" + "="*40)
    print("Case 1: c1 is not equal to 0")
    print("="*40)
    sympy.pprint(sympy.Eq(x_star, x_star_expr_c1_nonzero), use_unicode=True)

    print("\n" + "="*40)
    print("Case 2: c1 = 0")
    print("="*40)
    sympy.pprint(sympy.Eq(x_star, x_star_expr_c1_zero), use_unicode=True)

if __name__ == '__main__':
    find_transformation()