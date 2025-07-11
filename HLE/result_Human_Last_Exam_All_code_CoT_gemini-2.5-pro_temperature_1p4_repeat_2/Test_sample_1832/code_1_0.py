import sympy

def solve_christoffel_symbols():
    """
    Calculates and counts the non-zero Christoffel symbols for the Schwarzschild metric.
    This requires the sympy library. If you don't have it, please install it:
    pip install sympy
    """
    # Define coordinates and parameters
    t, r, theta, phi = sympy.symbols('t r theta phi')
    M = sympy.symbols('M', positive=True)

    # In geometric units, G=c=1. The Schwarzschild metric components are:
    # g_tt = -(1 - 2M/r)
    # g_rr = 1 / (1 - 2M/r)
    # g_thetatheta = r**2
    # g_phiphi = r**2 * sin(theta)**2

    # Metric tensor g_μν
    g = sympy.zeros(4, 4)
    g[0, 0] = -(1 - 2 * M / r)
    g[1, 1] = 1 / (1 - 2 * M / r)
    g[2, 2] = r**2
    g[3, 3] = r**2 * sympy.sin(theta)**2

    # Inverse metric tensor g^μν
    g_inv = g.inv()

    # List of coordinates
    coords = [t, r, theta, phi]
    coord_names = ['t', 'r', 'θ', 'φ']

    # Store non-zero symbols
    non_zero_symbols = []

    # Iterate over all unique components (μ <= ν)
    for rho in range(4):
        for mu in range(4):
            for nu in range(mu, 4):
                # Calculate Christoffel symbol Γ^ρ_μν
                gamma = sympy.simplify(sum([
                    sympy.Rational(1, 2) * g_inv[rho, sigma] * (
                        sympy.diff(g[nu, sigma], coords[mu]) +
                        sympy.diff(g[mu, sigma], coords[nu]) -
                        sympy.diff(g[mu, nu], coords[sigma])
                    )
                    for sigma in range(4)
                ]))

                if gamma != 0:
                    non_zero_symbols.append(((rho, mu, nu), gamma))
    
    print("The non-zero, unique Christoffel symbols Γ^ρ_{μν} (with μ <= ν) for the Schwarzschild metric are:")
    print("-" * 70)
    
    # Sort for consistent output
    non_zero_symbols.sort(key=lambda item: item[0])
    
    for (rho, mu, nu), value in non_zero_symbols:
        name = f"Γ^({coord_names[rho]})_({coord_names[mu]},{coord_names[nu]})"
        print(f"{name:<20} = {value}")

    print("-" * 70)

    # Count breakdown as requested
    counts_per_rho = {0: 0, 1: 0, 2: 0, 3: 0}
    for (rho, mu, nu), _ in non_zero_symbols:
        counts_per_rho[rho] += 1
    
    count_list = [counts_per_rho[i] for i in range(4)]
    total_count = len(non_zero_symbols)
    
    print("The number of unique non-zero Christoffel symbols is obtained by summing the counts for each upper index ρ:")
    final_equation = f"{count_list[0]} (for ρ=t) + {count_list[1]} (for ρ=r) + {count_list[2]} (for ρ=θ) + {count_list[3]} (for ρ=φ)"
    print(f"Total = {final_equation} = {total_count}")

# Execute the function
solve_christoffel_symbols()