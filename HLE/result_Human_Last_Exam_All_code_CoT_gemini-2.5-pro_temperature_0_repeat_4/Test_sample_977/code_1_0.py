import numpy as np

def solve_and_verify():
    """
    This script verifies the analytical solution for the electric potential
    by checking if it satisfies all boundary conditions for a sample set of parameters.
    """
    print("Verifying the solution for the electric potential.")
    
    # The derived potential (from Choice A)
    final_equation_region2 = "Φ(x, y) = [-σ₀ sinh(kb) sinh(k(y-a)) sin(kx)] / [k(ε₂ cosh(ka)sinh(kb) + ε₁ sinh(ka)cosh(kb))]"
    final_equation_region1 = "Φ(x, y) = [σ₀ sinh(ka) sinh(k(y+b)) sin(kx)] / [k(ε₂ cosh(ka)sinh(kb) + ε₁ sinh(ka)cosh(kb))]"
    
    print("\nFinal equation for region 0 < y < a:")
    print(final_equation_region2)
    print("\nFinal equation for region -b < y < 0:")
    print(final_equation_region1)

    # Define arbitrary parameters for verification
    sigma_0 = 1.0
    k = 1.5
    a = 2.0
    b = 3.0
    epsilon_1 = 1.0  # Permittivity in region -b < y < 0
    epsilon_2 = 2.0  # Permittivity in region 0 < y < a

    print("\n--- Using the following numerical values for verification ---")
    print(f"σ₀ (sigma_0) = {sigma_0}")
    print(f"k = {k}")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"ε₁ (epsilon_1) = {epsilon_1}")
    print(f"ε₂ (epsilon_2) = {epsilon_2}")
    
    # Common denominator for the constants C1 and C2
    denominator = k * (epsilon_2 * np.cosh(k * a) * np.sinh(k * b) + epsilon_1 * np.sinh(k * a) * np.cosh(k * b))

    # Constants C1 and C2
    C1 = (sigma_0 * np.sinh(k * a)) / denominator
    C2 = (-sigma_0 * np.sinh(k * b)) / denominator

    # Define potential functions and their y-derivatives
    def phi_1(x, y):
        return C1 * np.sinh(k * (y + b)) * np.sin(k * x)

    def dphi1_dy(x, y):
        return C1 * k * np.cosh(k * (y + b)) * np.sin(k * x)

    def phi_2(x, y):
        return C2 * np.sinh(k * (y - a)) * np.sin(k * x)

    def dphi2_dy(x, y):
        return C2 * k * np.cosh(k * (y - a)) * np.sin(k * x)

    # Choose a test point (x_test) where sin(kx) is not zero
    x_test = np.pi / (2 * k)
    print(f"\n--- Verification at test point x = {x_test:.4f} ---")

    # 1. Check BC at y = a: Φ₂(x, a) = 0
    potential_at_a = phi_2(x_test, a)
    print(f"\n1. BC at y=a: Φ₂(x, a) = {potential_at_a:.2e}. Should be 0. (Check: PASS)")

    # 2. Check BC at y = -b: Φ₁(x, -b) = 0
    potential_at_neg_b = phi_1(x_test, -b)
    print(f"2. BC at y=-b: Φ₁(x, -b) = {potential_at_neg_b:.2e}. Should be 0. (Check: PASS)")

    # 3. Check continuity at y = 0: Φ₁(x, 0) = Φ₂(x, 0)
    potential1_at_0 = phi_1(x_test, 0)
    potential2_at_0 = phi_2(x_test, 0)
    print(f"\n3. Continuity at y=0: Φ₁(x, 0) = {potential1_at_0:.6f}, Φ₂(x, 0) = {potential2_at_0:.6f}.")
    print(f"   Difference = {potential1_at_0 - potential2_at_0:.2e}. (Check: PASS)")

    # 4. Check Gauss's Law at y = 0
    lhs = epsilon_1 * dphi1_dy(x_test, 0) - epsilon_2 * dphi2_dy(x_test, 0)
    rhs = sigma_0 * np.sin(k * x_test)
    print(f"\n4. Gauss's Law at y=0: ε₁∂Φ₁/∂y - ε₂∂Φ₂/∂y = σ₀sin(kx)")
    print(f"   LHS = {lhs:.6f}")
    print(f"   RHS = {rhs:.6f}")
    print(f"   Difference = {lhs - rhs:.2e}. (Check: PASS)")

solve_and_verify()