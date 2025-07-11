import sympy

def solve_christoffel_symbols():
    """
    This function calculates the Christoffel symbols for the Schwarzschild metric
    and counts the number of non-zero components.
    """
    # 1. Define symbolic variables for coordinates and mass
    # Using t, r, th, ph for t, r, theta, phi
    t, r, th, ph, M = sympy.symbols('t r theta phi M')
    coords = [t, r, th, ph]

    # 2. Define the Schwarzschild metric tensor g_munu (with G=c=1)
    g = sympy.zeros(4)
    g[0, 0] = -(1 - 2 * M / r)
    g[1, 1] = 1 / (1 - 2 * M / r)
    g[2, 2] = r**2
    g[3, 3] = r**2 * sympy.sin(th)**2

    # 3. Calculate the inverse metric tensor g^munu
    g_inv = g.inv()
    # Simplify the inverse metric components for cleaner output
    g_inv.simplify()

    # 4. Calculate all Christoffel symbols
    # Initialize a 4x4x4 array to store the symbols
    christoffel = sympy.tensor.array.MutableDenseNDimArray([[[0]*4 for _ in range(4)] for _ in range(4)])
    
    # Pre-calculate derivatives of the metric tensor
    dg = sympy.tensor.array.derive_by_array(g, coords)

    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                s = 0
                for sigma in range(4):
                    # Christoffel symbol formula
                    term1 = dg[sigma, mu, nu]
                    term2 = dg[sigma, nu, mu]
                    term3 = dg[mu, nu, sigma]
                    s += g_inv[rho, sigma] * (term1 + term2 - term3)
                
                # Store the simplified result
                christoffel[rho, mu, nu] = sympy.simplify(s / 2)

    # 5. Count and print the non-zero symbols
    non_zero_count = 0
    print("The non-zero Christoffel symbols (Gamma^rho_mu,nu) and their values are:")
    print("Coordinates: 0=t, 1=r, 2=theta, 3=phi")
    print("-" * 60)
    
    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                if christoffel[rho, mu, nu] != 0:
                    non_zero_count += 1
                    # Print each non-zero symbol's indices and its value
                    print(f"Gamma^{rho}_{mu},{nu} = {christoffel[rho, mu, nu]}")

    # 6. Print the final count
    print("-" * 60)
    print(f"In total, there are {non_zero_count} non-zero Christoffel symbols.")

if __name__ == '__main__':
    solve_christoffel_symbols()