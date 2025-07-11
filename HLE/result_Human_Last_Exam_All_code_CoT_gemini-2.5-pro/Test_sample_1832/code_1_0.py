import sympy

def solve_christoffel_symbols():
    """
    Calculates the non-zero Christoffel symbols for the Schwarzschild metric.

    This function defines the Schwarzschild metric, calculates the Christoffel
    symbols using their formula, and then counts and prints the unique
    non-zero components.
    """
    # 1. Define symbolic variables for coordinates (t, r, theta, phi) and mass (M).
    # We use geometrized units where G = c = 1.
    t, r, theta, phi = sympy.symbols('t r theta phi')
    M = sympy.Symbol('M')
    coords = [t, r, theta, phi]
    coord_names = {t: 't', r: 'r', theta: '\\theta', phi: '\\phi'}

    # 2. Define the Schwarzschild metric tensor g_munu as a 4x4 matrix.
    # The metric is diagonal.
    g = sympy.zeros(4)
    schwarzschild_radius = 2 * M
    f = 1 - schwarzschild_radius / r
    g[0, 0] = -f
    g[1, 1] = 1/f
    g[2, 2] = r**2
    g[3, 3] = r**2 * sympy.sin(theta)**2

    # 3. Calculate the inverse metric tensor g_inv (g^rhosigma).
    g_inv = g.inv(method='ADJ') # Use adjoint method for symbolic matrices

    # 4. Calculate all first partial derivatives of the metric tensor.
    # dg[i][j,k] will store d(g_jk)/dx^i
    dg = [sympy.diff(g, coord) for coord in coords]

    # 5. Calculate all Christoffel symbols Gamma^rho_{mu,nu}.
    # The formula is Gamma^rho_{mu,nu} = 1/2 * g^rhosigma * (d_mu g_{sigma,nu} + d_nu g_{sigma,mu} - d_sigma g_{mu,nu})
    christoffel = sympy.MutableDenseNDimArray.zeros(4, 4, 4)
    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                sum_val = 0
                for sigma in range(4):
                    term1 = dg[mu][sigma, nu]  # d(g_sigmanu)/dx^mu
                    term2 = dg[nu][sigma, mu]  # d(g_sigmamu)/dx^nu
                    term3 = dg[sigma][mu, nu]  # d(g_munu)/dx^sigma
                    sum_val += 0.5 * g_inv[rho, sigma] * (term1 + term2 - term3)
                christoffel[rho, mu, nu] = sympy.simplify(sum_val)

    # 6. Count and display the unique non-zero symbols.
    # We account for the symmetry Gamma^rho_{mu,nu} = Gamma^rho_{nu,mu}.
    non_zero_symbols = {}
    for rho in range(4):
        for mu in range(4):
            for nu in range(mu, 4):  # Iterate mu <= nu to handle symmetry
                if christoffel[rho, mu, nu] != 0:
                    # Create a display-friendly name for the symbol
                    rho_name = coord_names[coords[rho]]
                    mu_name = coord_names[coords[mu]]
                    nu_name = coord_names[coords[nu]]
                    symbol_name = f"\\Gamma^{{{rho_name}}}_{{{mu_name}{nu_name}}}"
                    non_zero_symbols[symbol_name] = christoffel[rho, mu, nu]
    
    # 7. Print the results.
    print("The non-zero independent Christoffel symbols for the Schwarzschild metric are:")
    # The following equations show each non-zero symbol and its value.
    for name, value in sorted(non_zero_symbols.items()):
        # The 'equation' consists of the symbol name and its value.
        # We print both parts.
        print(f"{name} = {value}")

    print(f"\nIn total, there are {len(non_zero_symbols)} non-zero independent Christoffel symbols.")

if __name__ == '__main__':
    solve_christoffel_symbols()