from fractions import Fraction

def solve_policy_gradient():
    """
    Computes the policy gradients for the given MDP.
    """
    # --- Parameters Definition ---
    # Discount factor
    gamma = Fraction(1, 2)
    # Initial state distribution
    mu0_A = Fraction(1, 4)
    mu0_B = Fraction(3, 4)
    # Rewards r(s -> s')
    R_AA, R_AB = Fraction(1), Fraction(0)
    R_BA, R_BB = Fraction(1), Fraction(0)
    # Policy pi_1(s'|s)
    pi_A_to_A, pi_A_to_B = Fraction(1, 3), Fraction(2, 3)
    pi_B_to_A, pi_B_to_B = Fraction(1, 2), Fraction(1, 2)

    # --- Step 1: Calculate State-Value Functions V(s) ---
    # We solve the system of linear equations derived from the Bellman expectation equation:
    # V_s = sum_{s'} pi(s'|s) * [R(s,s') + gamma * V_{s'}]
    #
    # For state A: V_A = (1/3)*(1 + 1/2*V_A) + (2/3)*(0 + 1/2*V_B) => (5/6)V_A - (1/3)V_B = 1/3
    # For state B: V_B = (1/2)*(1 + 1/2*V_A) + (1/2)*(0 + 1/2*V_B) => (-1/4)V_A + (3/4)V_B = 1/2
    
    # Coefficients for the system [c11, c12; c21, c22] * [V_A; V_B] = [b1; b2]
    c11 = 1 - gamma * pi_A_to_A
    c12 = -gamma * pi_A_to_B
    c21 = -gamma * pi_B_to_A
    c22 = 1 - gamma * pi_B_to_B
    b1 = pi_A_to_A * R_AA + pi_A_to_B * R_AB
    b2 = pi_B_to_A * R_BA + pi_B_to_B * R_BB

    # Solve using Cramer's rule for V_A and V_B
    determinant = c11 * c22 - c12 * c21
    V_A = (b1 * c22 - b2 * c12) / determinant
    V_B = (b2 * c11 - b1 * c21) / determinant

    # --- Step 2: Calculate Action-Value Functions Q(s, a) ---
    # Q(s, s') = R(s, s') + gamma * V(s')
    Q_A_to_A = R_AA + gamma * V_A
    Q_A_to_B = R_AB + gamma * V_B

    # --- Step 3: Calculate Discounted State Visitation Frequency d(s) ---
    # We solve d^T = mu0^T * (I - gamma * P_pi)^-1
    # The matrix (I - gamma * P_pi) is the same one used for solving V(s).
    # We need its inverse, M_inv.
    # M_inv = (1/determinant) * [[c22, -c12], [-c21, c11]]
    M_inv_11 = c22 / determinant
    M_inv_12 = -c12 / determinant
    M_inv_21 = -c21 / determinant
    M_inv_22 = c11 / determinant
    
    # d_A = mu0_A * M_inv_11 + mu0_B * M_inv_21
    d_A = mu0_A * M_inv_11 + mu0_B * M_inv_21

    # --- Step 4: Compute the Gradients ---
    # Gradient formula: ∂V/∂π(a|s) = d(s) * Q(s, a)
    grad_pi_A_to_A = d_A * Q_A_to_A
    grad_pi_A_to_B = d_A * Q_A_to_B
    
    # --- Print Calculation and Results ---
    print("This problem requires computing policy gradients using the formula: ∂V(μ₀)/∂π(a|s) = d(s) * Q(s, a)\n")
    
    print("Step 1: Compute the value functions V(A) and V(B).")
    print(f"Solving the Bellman equations for the policy π₁ yields V(A) = {V_A.numerator}/{V_A.denominator} and V(B) = {V_B.numerator}/{V_B.denominator}.\n")
    
    print("Step 2: Compute the required Q-values for actions from state A.")
    print(f"Q(A, A->A) = r(A->A) + γ * V(A) = {R_AA} + ({gamma.numerator}/{gamma.denominator}) * ({V_A.numerator}/{V_A.denominator}) = {Q_A_to_A.numerator}/{Q_A_to_A.denominator}")
    print(f"Q(A, A->B) = r(A->B) + γ * V(B) = {R_AB} + ({gamma.numerator}/{gamma.denominator}) * ({V_B.numerator}/{V_B.denominator}) = {Q_A_to_B.numerator}/{Q_A_to_B.denominator}\n")
    
    print("Step 3: Compute the discounted state visitation frequency d(A).")
    print(f"By solving dᵀ = μ₀ᵀ(I - γP_π)⁻¹, we find d(A) = {d_A.numerator}/{d_A.denominator}.\n")
    
    print("Step 4: Compute the gradients.")
    print("Gradient 1: ∂V/∂π(A->A|A)")
    print(f"  = d(A) * Q(A, A->A)")
    print(f"  = ({d_A.numerator}/{d_A.denominator}) * ({Q_A_to_A.numerator}/{Q_A_to_A.denominator})")
    print(f"  = {grad_pi_A_to_A.numerator}/{grad_pi_A_to_A.denominator}\n")
    
    print("Gradient 2: ∂V/∂π(A->B|A)")
    print(f"  = d(A) * Q(A, A->B)")
    print(f"  = ({d_A.numerator}/{d_A.denominator}) * ({Q_A_to_B.numerator}/{Q_A_to_B.denominator})")
    print(f"  = {grad_pi_A_to_B.numerator}/{grad_pi_A_to_B.denominator}")

    # Final answer in the required format
    grad1_str = f"{grad_pi_A_to_A.numerator}/{grad_pi_A_to_A.denominator}"
    grad2_str = f"{grad_pi_A_to_B.numerator}/{grad_pi_A_to_B.denominator}"
    print(f"\n<<<{grad1_str}, {grad2_str}>>>")

solve_policy_gradient()