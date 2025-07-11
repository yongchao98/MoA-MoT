import sympy

def solve_christoffel_symbols():
    """
    Calculates and counts the non-zero Christoffel symbols for the Schwarzschild metric.
    """
    # Step 1: Define symbolic variables and coordinates
    # We use (t, r, theta, phi) as coordinates (0, 1, 2, 3)
    t, r, theta, phi = sympy.symbols('t r theta phi')
    M = sympy.Symbol('M')  # Mass of the gravitating body
    coords = [t, r, theta, phi]
    coord_names = ['t', 'r', 'theta', 'phi']

    # Step 2: Define the Schwarzschild metric tensor g_munu
    metric = sympy.zeros(4, 4)
    f = 1 - 2 * M / r
    metric[0, 0] = -f
    metric[1, 1] = 1 / f
    metric[2, 2] = r**2
    metric[3, 3] = r**2 * sympy.sin(theta)**2

    # Step 3: Calculate the inverse metric tensor g^rhosigma
    inverse_metric = metric.inv()

    # Step 4: Calculate the derivatives of the metric tensor
    # dg[k, mu, nu] = derivative of g_munu with respect to x^k
    dg = sympy.MutableDenseNDArray.zeros(4, 4, 4)
    for k in range(4):
        for mu in range(4):
            for nu in range(4):
                dg[k, mu, nu] = sympy.diff(metric[mu, nu], coords[k])

    # Step 5: Calculate all 64 Christoffel symbols Gamma^rho_munu
    christoffels = sympy.MutableDenseNDArray.zeros(4, 4, 4)
    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                sum_val = 0
                for sigma in range(4):
                    term = (dg[nu, sigma, mu] +
                            dg[mu, sigma, nu] -
                            dg[sigma, mu, nu])
                    sum_val += inverse_metric[rho, sigma] * term
                christoffels[rho, mu, nu] = sympy.simplify(0.5 * sum_val)

    # Step 6: Count the non-zero symbols, grouping by the upper index rho
    counts_per_rho = [0, 0, 0, 0]
    for rho in range(4):
        count = 0
        for mu in range(4):
            for nu in range(4):
                if christoffels[rho, mu, nu] != 0:
                    count += 1
        counts_per_rho[rho] = count
    
    # Step 7: Print the final result as an equation
    print("For a spherically symmetric gravitating body, the number of non-zero Christoffel symbols,")
    print(f"grouped by the upper index rho (for coordinates {coord_names}), are:")
    print(counts_per_rho)
    print("\nThe total number is the sum of these counts:")
    
    equation_parts = [str(c) for c in counts_per_rho]
    total_sum = sum(counts_per_rho)
    
    print(" + ".join(equation_parts) + f" = {total_sum}")

# Execute the function
solve_christoffel_symbols()