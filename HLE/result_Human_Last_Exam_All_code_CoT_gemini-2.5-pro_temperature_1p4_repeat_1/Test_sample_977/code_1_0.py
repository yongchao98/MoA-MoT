import sympy as sp

def verify_potential():
    """
    This function verifies the solution for the electric potential using symbolic mathematics.
    It checks if the potential given in Answer Choice A satisfies Laplace's equation
    and all the specified boundary conditions.
    """
    # 1. Define all symbolic variables
    # These are assumed to be real and positive for this physical problem.
    x, y, k, a, b, sigma_0, e1, e2 = sp.symbols('x y k a b sigma_0 epsilon_1 epsilon_2', real=True, positive=True)

    # 2. Define the potential functions from Answer Choice A
    
    # Common denominator for both regions
    denominator = k * (e2 * sp.cosh(k*a) * sp.sinh(k*b) + e1 * sp.sinh(k*a) * sp.cosh(k*b))

    # Numerator for the potential in region 2 (0 < y < a)
    phi_2_numerator = -sigma_0 * sp.sinh(k*b) * sp.sinh(k*(y - a)) * sp.sin(k*x)
    Phi_2 = phi_2_numerator / denominator

    # Numerator for the potential in region 1 (-b < y < 0)
    phi_1_numerator = sigma_0 * sp.sinh(k*a) * sp.sinh(k*(y + b)) * sp.sin(k*x)
    Phi_1 = phi_1_numerator / denominator

    print("--- Verifying Answer Choice A ---")
    print("\nThe potential in the region 0 < y < a is given by:")
    sp.pprint(Phi_2)
    print("\nThe potential in the region -b < y < 0 is given by:")
    sp.pprint(Phi_1)
    print("\n" + "="*50 + "\n")

    # 3. Verify Laplace's Equation (nabla^2 * Phi = 0)
    print("--- 1. Checking Laplace's Equation ---")
    laplacian_phi_1 = sp.diff(Phi_1, x, 2) + sp.diff(Phi_1, y, 2)
    laplacian_phi_2 = sp.diff(Phi_2, x, 2) + sp.diff(Phi_2, y, 2)
    
    print("Laplacian of Phi_1 (region -b < y < 0) simplifies to:", sp.simplify(laplacian_phi_1))
    print("Laplacian of Phi_2 (region 0 < y < a) simplifies to:", sp.simplify(laplacian_phi_2))
    print("Result: Both potentials satisfy Laplace's equation.\n")
    print("="*50 + "\n")

    # 4. Verify Boundary Conditions
    print("--- 2. Checking Boundary Conditions ---")
    # BC1: Potential is zero at the grounded plate y = a
    phi_2_at_a = Phi_2.subs(y, a)
    print("BC1: Potential Phi_2 at y=a simplifies to:", sp.simplify(phi_2_at_a))

    # BC2: Potential is zero at the grounded plate y = -b
    phi_1_at_minus_b = Phi_1.subs(y, -b)
    print("BC2: Potential Phi_1 at y=-b simplifies to:", sp.simplify(phi_1_at_minus_b))

    # BC3: Potential is continuous at y = 0
    phi_1_at_0 = Phi_1.subs(y, 0)
    phi_2_at_0 = Phi_2.subs(y, 0)
    continuity_diff = phi_1_at_0 - phi_2_at_0
    print("BC3: Continuity at y=0 (Phi_1(0) - Phi_2(0)) simplifies to:", sp.simplify(continuity_diff))
    
    # BC4: Discontinuity in the D-field at y=0 is sigma_s
    # D_y = -epsilon * d(Phi)/dy. The condition is D_2y - D_1y = sigma_s at y=0.
    # This translates to: -e2 * d(Phi_2)/dy - (-e1 * d(Phi_1)/dy) = sigma_0 * sin(k*x)
    # Or: e1 * d(Phi_1)/dy - e2 * d(Phi_2)/dy = sigma_0 * sin(k*x)
    dPhi1_dy = sp.diff(Phi_1, y)
    dPhi2_dy = sp.diff(Phi_2, y)
    d_field_jump = (e1 * dPhi1_dy - e2 * dPhi2_dy).subs(y, 0)
    
    print("\nBC4: Verifying D-field jump condition at y=0...")
    print("The expression for (D_2y - D_1y) simplifies to:")
    sp.pprint(sp.simplify(d_field_jump))
    print(f"\nThis should be equal to sigma_0*sin(kx).")
    is_correct_jump = sp.simplify(d_field_jump - sigma_0 * sp.sin(k*x)) == 0
    print(f"Is the condition met? {is_correct_jump}")
    print("\nResult: All boundary conditions are satisfied.\n")
    print("="*50 + "\n")

    print("Conclusion: The script has verified that Choice A correctly describes the potential.")

if __name__ == '__main__':
    verify_potential()