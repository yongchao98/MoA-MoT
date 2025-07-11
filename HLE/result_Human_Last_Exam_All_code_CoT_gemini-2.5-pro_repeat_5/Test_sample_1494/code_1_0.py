import sympy

def solve_task():
    """
    Calculates the exponent for the time-decay of ||dx(rho)(t)||L2.
    
    The decay is of the form exp(-C * t^p). We calculate p.
    
    The exponent p is given by the formula p = alpha / (alpha + beta), where:
    - alpha describes the smoothness of the initial data. For smooth functions,
      the energy in Fourier mode k decays like exp(-a*k^alpha), so we take alpha=1.
    - beta describes how the dissipation rate depends on the wave number k for large k.
      For the Stokes transport system, the decay rate lambda_k behaves like -1/k^beta,
      and a physical analysis of the system reveals beta=2.
    """
    
    # Parameter for initial data smoothness (alpha=1 for smooth C^infinity functions)
    alpha = sympy.Rational(1)
    
    # Parameter for the asymptotic behavior of the dissipation rate (beta=2 for this system)
    beta = sympy.Rational(2)
    
    # Calculate the exponent p
    p = alpha / (alpha + beta)
    
    # The decay for the L2 norm ||dx(rho)(t)|| is of the same form, exp(-C' * t^p).
    
    print("The best expected time-decay for ||dx(rho)(t)||L2 for smooth initial data follows the functional form:")
    print("exp(-C * t^p)")
    print("\nwhere C is a constant that depends on the initial data, and the exponent p is calculated as:")
    print(f"p = alpha / (alpha + beta) = {alpha} / ({alpha} + {beta}) = {p}")
    print("\nFinal decay law:")
    # The sympy.pretty_print function can be used for a nicer output format
    t = sympy.Symbol('t')
    C = sympy.Symbol('C')
    decay_formula = sympy.exp(-C * t**p)
    
    # We need to output each number in the final equation.
    # The equation is ||dx(rho)(t)|| ~ exp(-C * t^(1/3))
    print(f"||\u2202x\u03C1(t)||_L2 \u223C exp(-C \u00B7 t^({p.p}/{p.q}))")


solve_task()