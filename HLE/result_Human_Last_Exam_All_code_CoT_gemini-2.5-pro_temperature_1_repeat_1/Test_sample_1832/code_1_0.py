import sympy as sp

def calculate_christoffel_symbols():
    """
    This function symbolically calculates the non-zero Christoffel symbols
    for the Schwarzschild metric and prints the results.
    """
    # 1. Define symbolic variables and coordinates
    # We use t, r, theta, phi for coordinates and M for the mass of the body.
    # We work in units where G = c = 1.
    t, r, theta, phi = sp.symbols('t r theta phi')
    M = sp.Symbol('M')
    coords = [t, r, theta, phi]
    n = len(coords)

    # 2. Define the Schwarzschild metric tensor g_munu
    # The Schwarzschild radius is Rs = 2M.
    Rs = 2 * M
    g = sp.zeros(n)
    g[0, 0] = -(1 - Rs / r)
    g[1, 1] = 1 / (1 - Rs / r)
    g[2, 2] = r**2
    g[3, 3] = r**2 * sp.sin(theta)**2

    # 3. Calculate the inverse metric tensor g^munu
    g_inv = g.inv()

    # 4. Calculate all first partial derivatives of the metric tensor
    # dg[i][j][k] will store the derivative of g_ij with respect to the k-th coordinate.
    dg = [[[sp.diff(g[i, j], k) for k in coords] for j in range(n)] for i in range(n)]

    # 5. Calculate Christoffel symbols and count the non-zero ones
    print("The non-zero Christoffel symbols Gamma^rho_{mu,nu} are:")
    print("--------------------------------------------------------")
    
    counts_per_rho = [0] * n
    
    # Loop over all indices rho, mu, nu.
    # We use mu <= nu because the Christoffel symbol is symmetric in its lower indices.
    for rho in range(n):
        for mu in range(n):
            for nu in range(mu, n):
                gamma = sp.S(0)
                # Sum over the sigma index according to the formula
                for sigma in range(n):
                    term = g_inv[rho, sigma] * (dg[sigma][mu][nu] + dg[sigma][nu][mu] - dg[mu][nu][sigma])
                    gamma += term
                
                # The formula has a 1/2 factor
                gamma = sp.simplify(gamma / 2)

                # If the symbol is non-zero, print it and count it.
                if gamma != 0:
                    counts_per_rho[rho] += 1
                    rho_sym = coords[rho]
                    mu_sym = coords[mu]
                    nu_sym = coords[nu]
                    
                    print(f"Gamma^{rho_sym}_({mu_sym},{nu_sym}) = {gamma}")

    # 6. Print the final summary and the total count as an equation.
    print("\nSummary of unique non-zero components:")
    print("---------------------------------------")
    
    contributions = []
    for i in range(n):
        if counts_per_rho[i] > 0:
            contributions.append(str(counts_per_rho[i]))

    total_count = sum(counts_per_rho)
    equation = " + ".join(contributions)
    
    print("Number of non-zero symbols for each upper index (rho = t, r, theta, phi):")
    print(f"{counts_per_rho[0]} (for rho=t) + {counts_per_rho[1]} (for rho=r) + {counts_per_rho[2]} (for rho=theta) + {counts_per_rho[3]} (for rho=phi)")
    
    print("\nThe final equation for the total count is:")
    print(f"{equation} = {total_count}")

if __name__ == '__main__':
    calculate_christoffel_symbols()
    print("\n<<<9>>>")