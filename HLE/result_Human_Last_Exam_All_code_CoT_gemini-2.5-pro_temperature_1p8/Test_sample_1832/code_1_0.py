import sympy

def count_schwarzschild_christoffel_symbols():
    """
    This function calculates and counts the non-zero Christoffel symbols
    for the Schwarzschild metric.
    """
    # Define symbolic coordinates and constants
    # Using 'th' for theta and 'ph' for phi to avoid conflicts
    t, r, th, ph = sympy.symbols('t r theta phi')
    G, M = sympy.symbols('G M', positive=True)
    coords = [t, r, th, ph]
    coord_names = ['t', 'r', 'theta', 'phi']

    # Define the function f(r) for the Schwarzschild metric
    f = 1 - 2 * G * M / r

    # Define the Schwarzschild metric tensor g_munu (a 4x4 matrix)
    g_munu = sympy.Matrix([
        [-f, 0, 0, 0],
        [0, 1/f, 0, 0],
        [0, 0, r**2, 0],
        [0, 0, 0, r**2 * sympy.sin(th)**2]
    ])

    # Calculate the inverse metric tensor g^munu
    g_inv = g_munu.inv()

    # Calculate the partial derivatives of the metric tensor
    # g_diff[k][i][j] will store dg_{ij}/dx^k
    g_diff = [[[0 for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for k in range(4):
        for i in range(4):
            for j in range(4):
                g_diff[k][i][j] = sympy.diff(g_munu[i, j], coords[k])

    # Dictionary to store non-zero Christoffel symbols to avoid duplicates
    non_zero_symbols = {}

    # Calculate the Christoffel symbols Gamma^rho_mu_nu
    # The formula is: 1/2 * g^rho_sigma * (d_nu(g_sigma_mu) + d_mu(g_sigma_nu) - d_sigma(g_mu_nu))
    for rho in range(4):
        for mu in range(4):
            # Due to symmetry in lower indices (mu, nu), we only need to calculate for mu <= nu
            for nu in range(mu, 4):
                term = 0
                for sigma in range(4):
                    term += g_inv[rho, sigma] * (g_diff[nu][sigma][mu] + g_diff[mu][sigma][nu] - g_diff[sigma][mu][nu])
                
                # Simplify the final expression
                christoffel_symbol = sympy.simplify(0.5 * term)

                if christoffel_symbol != 0:
                    symbol_str = f"Γ^{coord_names[rho]}_{coord_names[mu]}{coord_names[nu]}"
                    equation = f"{symbol_str} = {christoffel_symbol}"
                    non_zero_symbols[symbol_str] = equation

    # Print the results
    print("The unique non-zero Christoffel symbols and their values are:\n")
    for symbol, equation in sorted(non_zero_symbols.items()):
        print(equation)
    
    count = len(non_zero_symbols)
    print(f"\nIn total, there are {count} unique non-zero Christoffel symbols for this metric.")

if __name__ == '__main__':
    count_schwarzschild_christoffel_symbols()