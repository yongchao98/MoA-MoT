import sympy

def solve_fluid_equation():
    """
    This function derives and displays the coefficients A(r) and B(r) for the
    governing differential equation of the fluid interface.

    The derivation is based on a variational approach (minimizing the total energy)
    and assumes that the surface tension gamma is a function of the radial
    position r, i.e., gamma(r), due to the influence of the electric field.
    """
    # Define symbols for the radial position 'r' and the surface tension 'gamma',
    # where gamma is a function of r.
    r = sympy.Symbol('r')
    gamma = sympy.Function('gamma')(r)

    # From the derivation based on the Euler-Lagrange equation for the system's
    # energy, the coefficient of the second derivative term is the surface tension itself.
    A_r = gamma

    # The coefficient of the first derivative term includes two effects:
    # 1. The geometry of the cylindrical system (gamma/r term from the Laplacian).
    # 2. The gradient of the surface tension (d(gamma)/dr term, a Marangoni-like effect).
    B_r = gamma / r + sympy.diff(gamma, r)

    # --- Output the results ---
    print("The governing linear differential equation for the interfacial shape xi(r) is of the form:")
    print("A(r) * d^2(xi)/dr^2 + B(r) * d(xi)/dr + C(r, xi) = 0")
    print("-" * 70)
    
    print("\nBased on the derivation that accounts for the electric field's influence on surface tension, the coefficients are:")
    
    # Print A(r)
    print("\nCoefficient A(r):")
    equation_A = sympy.Eq(sympy.Function('A')(r), A_r)
    sympy.pprint(equation_A, use_unicode=True)
    
    # Print B(r)
    print("\nCoefficient B(r):")
    equation_B = sympy.Eq(sympy.Function('B')(r), B_r)
    sympy.pprint(equation_B, use_unicode=True)


if __name__ == '__main__':
    solve_fluid_equation()