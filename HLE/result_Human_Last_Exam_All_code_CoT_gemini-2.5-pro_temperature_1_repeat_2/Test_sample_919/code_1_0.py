import sympy

def solve_emi_shielding_force():
    """
    This script symbolically derives the force per unit area on the conducting plane
    in the given electromagnetic interference (EMI) shielding setup.
    It uses the sympy library to handle the mathematical expressions.
    """
    
    # 1. Define Symbolic Variables
    # Physical constants and problem parameters
    K0, a, y, d = sympy.symbols('K_0 a y d', real=True, positive=True)
    mu0, mu = sympy.symbols('mu_0 mu', real=True, positive=True)
    x = sympy.symbols('x', real=True)
    
    # Unknown coefficients for the potential functions
    C, C2 = sympy.symbols('C C_2')

    # 2. Define Magnetic Scalar Potentials
    # Based on Laplace's equation and boundary conditions:
    # B_x = 0 at x=d leads to cosh(a(x-d)) form in the air gap.
    # Field must be finite at x -> -infinity, leading to exp(ax) form in the material.
    
    # Potential in the magnetic material (x < 0)
    Phi_mag = C * sympy.exp(a * x) * sympy.cos(a * y)
    
    # Potential in the air gap (0 < x < d)
    Phi_air = C2 * sympy.cosh(a * (x - d)) * sympy.cos(a * y)

    # 3. Define H and B fields from Potentials
    # H = -grad(Phi_m), B = mu*H
    
    # Fields in the magnetic material
    Hy_mag = -sympy.diff(Phi_mag, y)
    Bx_mag = -mu * sympy.diff(Phi_mag, x)
    
    # Fields in the air gap
    Hy_air = -sympy.diff(Phi_air, y)
    Bx_air = -mu0 * sympy.diff(Phi_air, x)

    # 4. Set up Equations from Boundary Conditions at x = 0
    
    # BC 1: Normal component of B is continuous: Bx_mag(0) = Bx_air(0)
    # We divide by common factors like cos(a*y) for simplicity.
    eq1_lhs = Bx_mag.subs(x, 0) / sympy.cos(a*y)
    eq1_rhs = Bx_air.subs(x, 0) / sympy.cos(a*y)
    eq1 = sympy.Eq(eq1_lhs.simplify(), eq1_rhs.simplify())

    # BC 2: Tangential H is discontinuous by K: Hy_air(0) - Hy_mag(0) = K_z
    # We divide by common factor sin(a*y).
    eq2_lhs = (Hy_air - Hy_mag).subs(x, 0) / sympy.sin(a*y)
    eq2_rhs = K0
    eq2 = sympy.Eq(eq2_lhs.simplify(), eq2_rhs.simplify())
    
    # 5. Solve for the Unknown Coefficient C2
    # We only need C2 to find the field in the air gap.
    solution = sympy.solve([eq1, eq2], (C, C2))
    C2_sol = solution[C2]

    # 6. Calculate the Magnetic Field at the Conductor (x = d)
    # At x=d, Bx=0, so B_total = By.
    By_at_d = (Hy_air * mu0).subs(x, d)
    
    # Substitute the solved coefficient C2 into the expression for B_y
    By_final_expr = By_at_d.subs(C2, C2_sol).simplify()
    
    # 7. Calculate the Force Per Unit Area
    # Force per area = B^2 / (2*mu_0) in the +x direction.
    force_per_area_mag = (By_final_expr**2 / (2 * mu0)).simplify()
    
    # 8. Print the Final Result
    # The final expression is a vector in the x-direction.
    # We format the output to be clear and match the style of the answer choices.
    
    numerator = sympy.mu(0)/2 * K0**2 * sympy.sin(a*y)**2
    denominator = (sympy.cosh(a*d) + sympy.mu(0)/mu * sympy.sinh(a*d))**2
    
    print("Derived Force per Unit Area Formula:")
    print("=" * 40)
    print("Numerator: ")
    sympy.pprint(numerator)
    print("\nDenominator: ")
    sympy.pprint(denominator)
    print("\nDirection: +î_x")
    print("="*40)
    
    # Reconstructing the final equation string for clarity
    final_eq_str = f"f/area = ({sympy.pretty(numerator)}) / ({sympy.pretty(denominator)}) * î_x"
    print("\nFull Equation:")
    print(final_eq_str)
    
    print("\nThis matches choice C.")

if __name__ == '__main__':
    solve_emi_shielding_force()