import sympy

def solve_lindhard_limit():
    """
    Calculates the Lindhard polarization function at zero frequency and zero momentum transfer.

    The static Lindhard function is Π(q, 0) = -g(E_F) * F(u), where g(E_F) is the density
    of states at the Fermi level and F(u) is a dimensionless function with u = q / (2*k_F).
    The value of Π(0, 0) is found by taking the limit of F(u) as u -> 0.
    """

    # Step 1: Define the dimensionless Lindhard function F(u) using sympy.
    # F(u) = 1/2 + (1 - u^2)/(4*u) * ln|(1 + u)/(1 - u)|
    u = sympy.Symbol('u')
    # We use sympy.log(1+u) - sympy.log(1-u) for ln((1+u)/(1-u))
    lindhard_func_dimensionless = 0.5 + (1 - u**2) / (4 * u) * (sympy.log(1 + u) - sympy.log(1 - u))

    # Step 2: Calculate the limit of the function as u approaches 0.
    # This corresponds to the zero momentum transfer (q -> 0) limit.
    limit_value = sympy.limit(lindhard_func_dimensionless, u, 0)

    # Step 3: Print the explanation and the final result.
    print("The Lindhard polarization function Π(q, ω) at zero frequency (ω=0) and zero momentum transfer (q=0) is sought.")
    print("The value is given by Π(0, 0) = -g(E_F) * F(0), where g(E_F) is the density of states at the Fermi level")
    print("and F(u) is the dimensionless Lindhard screening function, with u = q/(2*k_F).")
    print("We must compute the limit of F(u) as u -> 0.\n")
    print(f"The limit of the dimensionless part, lim_{{u->0}} F(u), is calculated to be: {limit_value}\n")

    print("This means the Lindhard function at zero momentum and frequency is equal to the negative")
    print("of the density of states at the Fermi level, as the dimensionless factor is 1.\n")

    # Final step: Output the final equation with each number explicitly shown.
    q_val = 0
    omega_val = 0
    # Convert from sympy's Integer type to a standard Python int for cleaner printing
    factor = int(limit_value)
    
    print("Final Equation:")
    print(f"Π(q={q_val}, ω={omega_val}) = -g(E_F) * {factor}")

solve_lindhard_limit()