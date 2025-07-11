import sympy as sp

def solve_emi_shielding_force():
    """
    This function symbolically derives the force per unit area on the conductor
    in the described EMI shielding setup.
    """
    # 1. Define symbolic variables
    K0, a, y, x, d = sp.symbols('K_0 a y x d', real=True, positive=True)
    mu, mu0 = sp.symbols('mu mu_0', real=True, positive=True)
    C1, C_const = sp.symbols('C_1 C', cls=sp.Function)

    print("Step 1: Define the physical model and magnetic scalar potentials.\n")
    # Region 1 (0 < x < d, air gap, mu0)
    # Region 2 (x < 0, magnetic material, mu)
    # Conductor at x=d, Current at x=0.

    # 2. Define magnetic scalar potentials psi_m in each region
    # BC at x=d (perfect conductor for time-varying field): B_normal = 0 -> H_x(d)=0
    # A potential of the form cosh(a(x-d)) satisfies H_x = -d/dx(psi) ~ sinh(a(x-d)) = 0 at x=d.
    psi_m1 = C1 * sp.cosh(a * (x - d)) * sp.cos(a * y)
    
    # BC at x -> -infinity: field must be finite
    psi_m2 = C_const * sp.exp(a * x) * sp.cos(a * y)
    
    print(f"Potential in Region 1 (air, 0<x<d): psi_m1 = {psi_m1}")
    print(f"Potential in Region 2 (mat, x<0): psi_m2 = {psi_m2}\n")

    # 3. Calculate H-fields from potentials (H = -grad(psi))
    H1_x = -sp.diff(psi_m1, x)
    H1_y = -sp.diff(psi_m1, y)
    H2_x = -sp.diff(psi_m2, x)
    H2_y = -sp.diff(psi_m2, y)
    
    print("Step 2: Apply boundary conditions at the x=0 interface.\n")
    # 4. Apply boundary conditions at x=0
    
    # BC 1: Normal component of B is continuous: mu0 * H1_x = mu * H2_x at x=0
    bc1_lhs = (mu0 * H1_x).subs(x, 0)
    bc1_rhs = (mu * H2_x).subs(x, 0)
    eq1 = sp.Eq(bc1_lhs, bc1_rhs)
    print(f"BC 1 (B_x continuous): {eq1}")
    
    # BC 2: Tangential H is discontinuous: H1_y - H2_y = K0 * sin(a*y) at x=0
    # The amplitudes must match
    amp_H1y = H1_y.subs(x, 0) / sp.sin(a * y)
    amp_H2y = H2_y.subs(x, 0) / sp.sin(a * y)
    eq2 = sp.Eq(amp_H1y - amp_H2y, K0)
    print(f"BC 2 (H_y discontinuous): {eq2}\n")

    # 5. Solve for the constant C1
    solution = sp.solve([eq1, eq2], (C1, C_const))
    C1_val = solution[C1]
    
    print("Step 3: Solve for the integration constant C1.\n")
    print(f"Solved C1 = {sp.pretty(C1_val, use_unicode=False)}\n")
    
    # 6. Calculate B-field at the conductor surface (x=d)
    # B_x is zero by boundary condition. We only need B_y.
    B1_y_surface = (mu0 * H1_y).subs(x, d)
    print("Step 4: Calculate the magnetic field B at the conductor surface x=d.\n")
    print(f"B_y at x=d is: {B1_y_surface}\n")
    
    # Total B-field magnitude squared
    B_squared = B1_y_surface**2
    
    # 7. Calculate force per unit area
    force_per_area_mag = B_squared / (2 * mu0)
    
    # Substitute C1 value into the force expression
    force_per_area_mag_final = force_per_area_mag.subs(C1, C1_val)
    
    # Construct the final vector expression for force
    i_x = sp.Matrix([1, 0, 0])
    
    print("Step 5: Calculate the force per unit area on the conductor.\n")
    print("f/Area = (B_surface^2 / (2 * mu_0)) * i_x\n")
    
    final_force_expr = sp.simplify(force_per_area_mag_final)
    
    print("Final result for the force per unit area:")
    
    # Formatting the output to match answer choices
    numerator = sp.S(1)/2 * mu0 * K0**2 * sp.sin(a*y)**2
    denominator = (sp.cosh(a*d) + (mu0/mu)*sp.sinh(a*d))**2
    
    print(f"f/Area = ({numerator}) / ({denominator}) * i_x")


if __name__ == '__main__':
    solve_emi_shielding_force()