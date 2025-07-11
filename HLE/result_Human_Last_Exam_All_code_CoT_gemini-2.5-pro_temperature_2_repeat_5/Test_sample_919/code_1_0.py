import sympy as sp

def solve_force_on_conductor():
    """
    Symbolically calculates the force per unit area on the conducting plane
    using the defined physical system and boundary conditions.
    """
    # Define the symbolic variables used in the problem
    x, y, d, a = sp.symbols('x y d a', real=True, positive=True)
    K0, mu, mu0 = sp.symbols('K_0 mu mu_0', real=True, positive=True)
    C1, C2 = sp.symbols('C1 C2', real=True)

    # Define the magnetic scalar potential (Φ_m) for each region
    # The form of the potential is derived from Laplace's equation and boundary conditions.
    
    # Region 1 (x < 0, magnetic material)
    # The term e^(-ax) is dropped so the field doesn't diverge at x -> -∞
    Phi1 = C1 * sp.exp(a * x) * sp.cos(a * y)
    
    # Region 2 (0 < x < d, air gap)
    # The form C2*cosh(a*(x-d)) is chosen to automatically satisfy the boundary
    # condition B_x = 0 at the perfect conductor surface x = d.
    Phi2 = C2 * sp.cosh(a * (x - d)) * sp.cos(a * y)

    # Derive the H-field components (H = -∇Φ_m)
    H1x = -sp.diff(Phi1, x)
    H1y = -sp.diff(Phi1, y)
    H2x = -sp.diff(Phi2, x)
    H2y = -sp.diff(Phi2, y)

    # Apply boundary conditions at the x = 0 interface to create a system of equations
    # BC 1: The normal component of B is continuous (μ*H_{1x} = μ₀*H_{2x})
    eq1 = sp.Eq(mu * H1x.subs(x, 0), mu0 * H2x.subs(x, 0))
    
    # BC 2: The tangential component of H is discontinuous by the surface current K_z
    # H_{2y} - H_{1y} = K_z = K_0*sin(ay)
    eq2 = sp.Eq(H2y.subs(x, 0) - H1y.subs(x, 0), K0 * sp.sin(a * y))

    # Solve the system of equations for the unknown coefficients C1 and C2
    # We remove the trivial trigonometric factors to simplify solving
    eq1_simplified = sp.Eq(eq1.lhs / sp.cos(a*y), eq1.rhs / sp.cos(a*y))
    eq2_simplified = sp.Eq(eq2.lhs / sp.sin(a*y), eq2.rhs / sp.sin(a*y))
    solution = sp.solve([eq1_simplified, eq2_simplified], (C1, C2))
    
    # The coefficient C2 is needed for the field in the air gap
    C2_solved = solution[C2]

    # Calculate the magnetic field B parallel to the surface at x = d
    # H_y(d) = -∂Φ_2/∂y at x=d
    Hy_at_d = H2y.subs(x, d).subs(C2, C2_solved)
    # B_y = μ₀ * H_y in the air gap
    By_at_d = mu0 * Hy_at_d

    # Calculate the force per unit area on the conductor.
    # The force is due to magnetic pressure P = B_{||}^2 / (2μ₀) and acts in the +x direction.
    force_per_area = By_at_d**2 / (2 * mu0)
    
    # Simplify the final expression
    force_per_area_simplified = sp.simplify(force_per_area)
    
    # The problem asks for the vector expression of the force. The force is in the x-direction.
    force_vector_str = f"({force_per_area_simplified}) * i_x"
    
    print("The symbolic expression for the force per unit y-z area is:")
    print(force_vector_str)
    
    # Let's reconstruct the result in a more readable format similar to the options
    final_numerator = mu0 * K0**2 * sp.sin(a*y)**2
    final_denominator = 2 * (sp.cosh(a*d) + (mu0/mu)*sp.sinh(a*d))**2
    pretty_expression = final_numerator / final_denominator
    
    print("\nIn formatted mathematical notation:")
    sp.pprint(pretty_expression)
    print(" "*25 + " * î_x")

if __name__ == '__main__':
    solve_force_on_conductor()
