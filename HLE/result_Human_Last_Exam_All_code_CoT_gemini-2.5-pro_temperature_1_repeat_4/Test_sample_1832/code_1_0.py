import sympy

def solve_christoffel_symbols():
    """
    Calculates and counts the non-zero Christoffel symbols for the Schwarzschild metric.
    """
    # Define coordinates and the mass parameter M. We use units where G=c=1.
    t, r, theta, phi = sympy.symbols('t r theta phi')
    M = sympy.Symbol('M')
    
    # Coordinates list
    coords = [t, r, theta, phi]
    coord_names = ['t', 'r', 'theta', 'phi']

    # Define the Schwarzschild metric components
    f_r = (1 - 2*M / r)
    g_tt = -f_r
    g_rr = 1 / f_r
    g_thetatheta = r**2
    g_phiphi = r**2 * sympy.sin(theta)**2

    # Create the metric tensor as a sympy Matrix
    g = sympy.Matrix([
        [g_tt, 0, 0, 0],
        [0, g_rr, 0, 0],
        [0, 0, g_thetatheta, 0],
        [0, 0, 0, g_phiphi]
    ])

    # Calculate the inverse metric tensor
    g_inv = g.inv()

    # Calculate the partial derivatives of the metric tensor with respect to all coordinates
    # dg[i, j, k] will store dg_{jk} / dx^i
    dg = sympy.tensor.array.derive_by_array(g, coords)

    # Initialize a 4x4x4 array for the Christoffel symbols
    christoffel = sympy.MutableDenseNDimArray.zeros(4, 4, 4)

    # Calculate Christoffel symbols using the formula
    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                term = sympy.S.Zero
                for sigma in range(4):
                    # Formula: 1/2 * g^{rho,sigma} * (d_mu g_{nu,sigma} + d_nu g_{mu,sigma} - d_sigma g_{mu,nu})
                    term += g_inv[rho, sigma] * (dg[mu, nu, sigma] + dg[nu, mu, sigma] - dg[sigma, mu, nu])
                
                christoffel[rho, mu, nu] = sympy.simplify(sympy.S.Half * term)

    # Count and print the unique non-zero Christoffel symbols
    print("The non-zero Christoffel symbols (and their corresponding equations) for the Schwarzschild metric are:\n")
    
    unique_non_zero_count = 0
    
    # Iterate over unique components (mu <= nu) due to symmetry
    for rho in range(4):
        for mu in range(4):
            for nu in range(mu, 4):
                symbol = christoffel[rho, mu, nu]
                if symbol != 0:
                    unique_non_zero_count += 1
                    
                    # Format the symbol's name for printing
                    rho_name = coord_names[rho]
                    mu_name = coord_names[mu]
                    nu_name = coord_names[nu]
                    
                    base_name = f"Gamma^{rho_name}_{mu_name}{nu_name}"
                    
                    # If mu != nu, show the symmetric component as well
                    if mu != nu:
                        symmetric_name = f"Gamma^{rho_name}_{nu_name}{mu_name}"
                        print(f"{base_name} = {symmetric_name} = {symbol}")
                    else:
                        print(f"{base_name} = {symbol}")

    print(f"\nTotal number of unique non-zero Christoffel symbols: {unique_non_zero_count}")

if __name__ == '__main__':
    solve_christoffel_symbols()