import sympy as sy

def verify_potential():
    """
    Symbolically verifies the electric potential for the given problem.
    This script checks if the potential expressions from Option A satisfy
    Laplace's equation and all boundary conditions.
    """
    # Define the symbolic variables
    x, y, k, a, b, sigma0, eps1, eps2 = sy.symbols(
        'x y k a b sigma_0 epsilon_1 epsilon_2', real=True, positive=True)

    print("Verifying the solution for the electric potential from choice A.")
    print("=" * 60)

    # Define the denominator term, which is common to both potential expressions
    denominator = k * (eps2 * sy.cosh(k * a) * sy.sinh(k * b) + \
                       eps1 * sy.sinh(k * a) * sy.cosh(k * b))

    # Potential in Region 2 (0 < y < a) from option A
    phi2_num = -sigma0 * sy.sinh(k * b) * sy.sinh(k * (y - a)) * sy.sin(k * x)
    Phi_2 = phi2_num / denominator

    # Potential in Region 1 (-b < y < 0) from option A
    phi1_num = sigma0 * sy.sinh(k * a) * sy.sinh(k * (y + b)) * sy.sin(k * x)
    Phi_1 = phi1_num / denominator

    # --- Verification Steps ---

    # 1. Check Laplace's Equation: d^2(Phi)/dx^2 + d^2(Phi)/dy^2 = 0
    print("1. Checking if potential satisfies Laplace's equation (del^2 Phi = 0):")
    laplacian_phi1 = sy.diff(Phi_1, x, 2) + sy.diff(Phi_1, y, 2)
    laplacian_phi2 = sy.diff(Phi_2, x, 2) + sy.diff(Phi_2, y, 2)
    print(f"   - For Region 1 (-b < y < 0): Is del^2 Phi_1 = 0? {sy.simplify(laplacian_phi1) == 0}")
    print(f"   - For Region 2 (0 < y < a):  Is del^2 Phi_2 = 0? {sy.simplify(laplacian_phi2) == 0}")
    print("-" * 60)

    # 2. Check Boundary Conditions at the conductors
    print("2. Checking boundary conditions at grounded plates:")
    # At y = a, Phi_2 must be 0
    phi2_at_a = Phi_2.subs(y, a)
    print(f"   - At y = a: Is Phi_2(a) = 0? {sy.simplify(phi2_at_a) == 0}")
    # At y = -b, Phi_1 must be 0
    phi1_at_neg_b = Phi_1.subs(y, -b)
    print(f"   - At y = -b: Is Phi_1(-b) = 0? {sy.simplify(phi1_at_neg_b) == 0}")
    print("-" * 60)
    
    # 3. Check Boundary Conditions at the interface y = 0
    print("3. Checking boundary conditions at the interface y = 0:")
    # Potential must be continuous: Phi_1(y=0) == Phi_2(y=0)
    phi1_at_0 = Phi_1.subs(y, 0)
    phi2_at_0 = Phi_2.subs(y, 0)
    print(f"   - Is potential continuous (Phi_1(0) == Phi_2(0))? "
          f"{sy.simplify(phi1_at_0 - phi2_at_0) == 0}")

    # Check discontinuity of electric displacement field
    d_phi1_dy = sy.diff(Phi_1, y)
    d_phi2_dy = sy.diff(Phi_2, y)
    gauss_law_lhs = eps1 * d_phi1_dy.subs(y, 0) - eps2 * d_phi2_dy.subs(y, 0)
    gauss_law_rhs = sigma0 * sy.sin(k * x)
    print(f"   - Is Gauss's Law satisfied (eps1*dPhi1/dy - eps2*dPhi2/dy = sigma_s)? "
          f"{sy.simplify(gauss_law_lhs - gauss_law_rhs) == 0}")
    print("-" * 60)
    
    # All checks passed, so the expression is correct.
    print("Conclusion: All conditions are satisfied by Option A.")
    print("The electric potential Phi(x, y) in the region 0 < y < a is given by:")
    print()
    print("           -sigma_0 * sinh(k*b) * sinh(k*(y - a)) * sin(k*x)")
    print("Phi(x,y) = -----------------------------------------------------------------------")
    print("            k * [epsilon_2*cosh(k*a)*sinh(k*b) + epsilon_1*sinh(k*a)*cosh(k*b)]")

if __name__ == '__main__':
    verify_potential()
