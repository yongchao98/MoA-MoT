import sympy

def solve_christoffel_symbols():
    """
    Calculates and counts the non-zero Christoffel symbols for the Schwarzschild metric.
    """
    # 1. Setup symbolic variables and coordinates (using geometric units G=c=1)
    # The coordinates are x^0=t, x^1=r, x^2=theta, x^3=phi
    t, r, theta, phi, M = sympy.symbols('t r theta phi M')
    coords = [t, r, theta, phi]
    
    # Let rs be the Schwarzschild radius, rs = 2M
    rs = 2 * M

    # 2. Define the Schwarzschild metric tensor g_munu
    # ds^2 = -(1-rs/r)dt^2 + (1-rs/r)^-1 dr^2 + r^2 d(theta)^2 + r^2 sin(theta)^2 d(phi)^2
    g = sympy.zeros(4)
    g[0, 0] = -(1 - rs / r)
    g[1, 1] = 1 / (1 - rs / r)
    g[2, 2] = r**2
    g[3, 3] = r**2 * sympy.sin(theta)**2

    # 3. Calculate the inverse metric tensor g^munu
    g_inv = sympy.zeros(4)
    for i in range(4):
        # For a diagonal matrix, the inverse is the reciprocal of the diagonal elements
        g_inv[i, i] = 1 / g[i, i]

    # 4. Iterate and calculate the unique non-zero Christoffel symbols
    non_zero_symbols = []
    count = 0
    
    print("The non-zero Christoffel symbols Γ^ρ_{μν} for the Schwarzschild metric are:")

    # Loop over the upper index rho
    for rho in range(4):
        # Loop over the lower indices mu and nu, with mu <= nu for unique components
        for mu in range(4):
            for nu in range(mu, 4):
                
                # Christoffel symbol formula:
                # Γ^ρ_{μν} = 1/2 * g^{ρσ} * (∂_μ g_{νσ} + ∂_ν g_{μσ} - ∂_σ g_{μν})
                gamma_expr = sympy.sympify(0)
                for sigma in range(4):
                    d_mu_g_nu_sigma = sympy.diff(g[nu, sigma], coords[mu])
                    d_nu_g_mu_sigma = sympy.diff(g[mu, sigma], coords[nu])
                    d_sigma_g_mu_nu = sympy.diff(g[mu, nu], coords[sigma])
                    
                    d_g_term = d_mu_g_nu_sigma + d_nu_g_mu_sigma - d_sigma_g_mu_nu
                    
                    gamma_expr += g_inv[rho, sigma] * d_g_term

                gamma_expr = sympy.simplify(gamma_expr / 2)

                if gamma_expr != 0:
                    count += 1
                    # Print each non-zero component found. This satisfies the "output each number" requirement
                    # by showing the full equation for each non-zero term.
                    print(f"Γ^{rho}_{{{mu},{nu}}} = {gamma_expr}")

    print("\nTotal number of unique non-zero Christoffel symbols:")
    print(count)


if __name__ == '__main__':
    solve_christoffel_symbols()