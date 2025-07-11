import sympy
from sympy import sin, cos, symbols, simplify, pprint

def solve_and_verify():
    """
    This script verifies the solution to the magnetostatics problem
    by checking the boundary conditions for the proposed fields.
    """
    # Define symbolic variables for the problem
    r, theta = symbols('r, theta', real=True, positive=True)
    M0, Rp, R = symbols('M_0 R_p R', real=True, positive=True)

    print("This program will verify the correct answer (Choice B) by checking the boundary conditions.")
    print("-" * 50)
    print("The proposed solution is:")
    
    # --- Define H fields and M vector from Answer Choice B ---

    # Region 1 (0 < r < Rp)
    H1_coeff = M0 * (2*Rp**3 + R**3) / (3*R**3)
    H1r_expr = H1_coeff * (-cos(theta))
    H1t_expr = H1_coeff * (sin(theta))
    
    print("\nIn the region 0 < r < R_p:")
    print("H_r = M_0 * (2*R_p**3 + R**3)/(3*R**3) * (-cos(theta))")
    print("H_theta = M_0 * (2*R_p**3 + R**3)/(3*R**3) * (sin(theta))")

    # Region 2 (Rp < r < R)
    term_R = (Rp/R)**3
    term_r = (Rp/r)**3
    H2r_expr = -2*M0/3 * (term_R - term_r) * cos(theta)
    H2t_expr = M0/3 * (2*term_R + term_r) * sin(theta)

    print("\nIn the region R_p < r < R:")
    print("H_r = -2*M_0/3 * [ (R_p/R)**3 - (R_p/r)**3 ] * cos(theta)")
    print("H_theta = M_0/3 * [ 2*(R_p/R)**3 + (R_p/r)**3 ] * sin(theta)")
    
    # Magnetization M vector in spherical components
    Mr_expr = M0 * cos(theta)
    Mt_expr = -M0 * sin(theta)

    # --- Define B fields (setting mu0 = 1 for simplicity in boundary checks) ---
    # B = mu0 * (H + M)
    B1r_expr = H1r_expr + Mr_expr
    B2r_expr = H2r_expr

    # --- Verification Step 1: Boundary condition at r = R (perfect conductor) ---
    print("\n" + "-"*50)
    print("1. VERIFICATION at r = R (Perfect Conductor)")
    print("   Condition: The normal component of the B-field (and H-field) must be zero.")
    
    H2r_at_R = H2r_expr.subs(r, R)
    is_bc1_satisfied = (simplify(H2r_at_R) == 0)
    print(f"   H_r in Region 2 evaluated at r=R simplifies to: {simplify(H2r_at_R)}")
    print(f"   --> Boundary condition at r=R satisfied: {is_bc1_satisfied}")

    # --- Verification Step 2: Boundary conditions at r = Rp (interface) ---
    print("\n" + "-"*50)
    print("2. VERIFICATION at r = R_p (Interface)")
    
    # a) Continuity of tangential H
    print("\n   a) Condition: The tangential component of H (H_theta) must be continuous.")
    H1t_at_Rp = H1t_expr.subs(r, Rp)
    H2t_at_Rp = H2t_expr.subs(r, Rp)
    diff_Ht = simplify(H1t_at_Rp - H2t_at_Rp)
    is_bc2a_satisfied = (diff_Ht == 0)

    print("      H_theta in Region 1 at r=R_p:", simplify(H1t_at_Rp))
    print("      H_theta in Region 2 at r=R_p:", simplify(H2t_at_Rp))
    print(f"      Difference (H1_theta - H2_theta) at r=R_p simplifies to: {diff_Ht}")
    print(f"      --> H_theta continuity satisfied: {is_bc2a_satisfied}")

    # b) Continuity of normal B
    print("\n   b) Condition: The normal component of B (B_r) must be continuous.")
    B1r_at_Rp = B1r_expr.subs(r, Rp)
    B2r_at_Rp = B2r_expr.subs(r, Rp)
    diff_Br = simplify(B1r_at_Rp - B2r_at_Rp)
    is_bc2b_satisfied = (diff_Br == 0)
    
    print("      B_r in Region 1 at r=R_p simplifies to:", simplify(B1r_at_Rp))
    print("      B_r in Region 2 at r=R_p simplifies to:", simplify(B2r_at_Rp))
    print(f"      Difference (B1_r - B2_r) at r=R_p simplifies to: {diff_Br}")
    print(f"      --> B_r continuity satisfied: {is_bc2b_satisfied}")
    
    print("\n" + "-"*50)
    if is_bc1_satisfied and is_bc2a_satisfied and is_bc2b_satisfied:
        print("CONCLUSION: All boundary conditions are satisfied. Choice B is the correct solution.")
    else:
        print("CONCLUSION: One or more boundary conditions failed. The solution is incorrect.")

# Execute the verification
solve_and_verify()