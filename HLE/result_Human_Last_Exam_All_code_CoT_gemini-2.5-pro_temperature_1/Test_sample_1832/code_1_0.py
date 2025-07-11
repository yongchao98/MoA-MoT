import sympy

def calculate_christoffel_symbols():
    """
    Calculates and prints the non-zero Christoffel symbols for the Schwarzschild metric.
    """
    # Define coordinates and parameters using sympy
    # We use th and ph for theta and phi to avoid conflicts with sympy functions
    coords = t, r, th, ph = sympy.symbols('t r \\theta \\phi')
    M = sympy.Symbol('M')

    # Schwarzschild radius in geometrized units (G=c=1)
    rs = 2 * M

    # Define the Schwarzschild metric tensor components
    g_tt = -(1 - rs / r)
    g_rr = 1 / (1 - rs / r)
    g_thth = r**2
    g_phph = r**2 * sympy.sin(th)**2

    # Create the metric tensor as a 4x4 matrix
    g = sympy.Matrix([
        [g_tt, 0, 0, 0],
        [0, g_rr, 0, 0],
        [0, 0, g_thth, 0],
        [0, 0, 0, g_phph]
    ])

    # Calculate the inverse of the metric tensor
    g_inv = g.inv()
    g_inv.simplify()

    # Pre-calculate all partial derivatives of the metric tensor components
    # dg[i][j][k] will store the derivative of g_{ij} with respect to the k-th coordinate
    dg = [[[sympy.diff(g[i, j], k) for k in coords] for j in range(4)] for i in range(4)]

    # Dictionary to store non-zero symbols, grouped by the upper index (rho)
    non_zero_symbols = {}

    # Iterate through all indices rho, mu, nu to calculate Gamma^rho_{mu,nu}
    for rho in range(4):
        non_zero_symbols[rho] = []
        for mu in range(4):
            # Due to symmetry Gamma^rho_{mu,nu} = Gamma^rho_{nu,mu}, we only calculate for mu <= nu
            for nu in range(mu, 4):
                gamma_sum = 0
                # Sum over the index sigma
                for sigma in range(4):
                    term1 = dg[nu][sigma][mu]  # d(g_{nu,sigma}) / dx^mu
                    term2 = dg[mu][sigma][nu]  # d(g_{mu,sigma}) / dx^nu
                    term3 = dg[mu][nu][sigma]  # d(g_{mu,nu}) / dx^sigma
                    gamma_sum += g_inv[rho, sigma] * (term1 + term2 - term3)

                # The formula has a factor of 1/2
                gamma = sympy.simplify(gamma_sum / 2)

                if gamma != 0:
                    # Store the non-zero symbol
                    non_zero_symbols[rho].append((mu, nu, gamma))

    # --- Output the results ---
    print("The non-zero Christoffel symbols are:")
    print("-" * 40)
    coord_names = ['t', 'r', '\\theta', '\\phi']
    total_count = 0

    for rho, symbols_list in non_zero_symbols.items():
        if not symbols_list:
            continue
        for mu, nu, expr in symbols_list:
            # Format the symbol name for printing
            upper_idx = coord_names[rho]
            lower_idx1 = coord_names[mu]
            lower_idx2 = coord_names[nu]
            print(f"Gamma^{upper_idx}_{{{lower_idx1}{lower_idx2}}} = {expr}")
        total_count += len(symbols_list)
    print("-" * 40)
    
    print("The final equation for the total number of unique non-zero symbols is:")
    count_per_rho = [len(non_zero_symbols.get(i, [])) for i in range(4)]
    equation_str = " + ".join(map(str, count_per_rho))
    print(f"Total Symbols = {equation_str} = {total_count}")


if __name__ == '__main__':
    calculate_christoffel_symbols()