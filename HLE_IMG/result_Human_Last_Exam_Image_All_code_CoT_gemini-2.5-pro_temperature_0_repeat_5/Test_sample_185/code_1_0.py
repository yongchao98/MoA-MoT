import numpy as np
from fractions import Fraction

def solve_policy_gradient():
    """
    Computes the policy gradients for the given MDP.
    """
    # --- 1. Define MDP and Policy Parameters ---
    gamma = Fraction(1, 2)

    # Initial distribution mu_0
    mu0 = np.array([Fraction(1, 4), Fraction(3, 4)]) # [mu0(A), mu0(B)]

    # Policy pi_1 probabilities
    p_AA = Fraction(1, 3)
    p_AB = Fraction(2, 3)
    p_BA = Fraction(1, 2)
    p_BB = Fraction(1, 2)

    # Rewards r(s, a, s')
    r_A_to_A = 1
    r_A_to_B = 0
    r_B_to_A = 1
    r_B_to_B = 0

    # --- 2. Calculate State-Value Function V^pi ---
    # The Bellman equation is V = R_pi + gamma * P_pi * V
    # which can be rewritten as (I - gamma * P_pi) * V = R_pi
    
    # P_pi is the transition matrix under the policy
    P_pi = np.array([
        [p_AA, p_AB],  # Transitions from A
        [p_BA, p_BB]   # Transitions from B
    ])

    # R_pi is the expected immediate reward from each state
    R_pi_A = p_AA * r_A_to_A + p_AB * r_A_to_B
    R_pi_B = p_BA * r_B_to_A + p_BB * r_B_to_B
    R_pi = np.array([R_pi_A, R_pi_B])

    # Solve the system for V
    I = np.identity(2)
    A_V = I - gamma * P_pi
    V_pi = np.linalg.solve(A_V.astype(float), R_pi.astype(float))
    V_A, V_B = V_pi[0], V_pi[1]

    # --- 3. Calculate Action-Value Function Q^pi ---
    # Q(s, a) = R(s, a) + gamma * V(s')
    Q_A_to_A = r_A_to_A + gamma * V_A
    Q_A_to_B = r_A_to_B + gamma * V_B

    # --- 4. Calculate Discounted State Visitation Frequencies d ---
    # The equation for d is (I - gamma * P_pi)^T * d = mu0
    A_d = (I - gamma * P_pi).T
    d = np.linalg.solve(A_d.astype(float), mu0.astype(float))
    d_A, d_B = d[0], d[1]

    # --- 5. Compute the Gradients ---
    grad_p_AA = d_A * Q_A_to_A
    grad_p_AB = d_A * Q_A_to_B

    # --- 6. Print the results ---
    # Convert floats to fractions for clear, exact representation
    f_d_A = Fraction(d_A).limit_denominator()
    f_Q_A_to_A = Fraction(Q_A_to_A).limit_denominator()
    f_grad_p_AA = Fraction(grad_p_AA).limit_denominator()
    
    f_Q_A_to_B = Fraction(Q_A_to_B).limit_denominator()
    f_grad_p_AB = Fraction(grad_p_AB).limit_denominator()

    print("Gradient of V w.r.t. pi(A->A|A):")
    print(f"dV/d(pi(A->A|A)) = d_mu(A) * Q_pi(A, A->A)")
    print(f"                 = {f_d_A} * {f_Q_A_to_A}")
    print(f"                 = {f_grad_p_AA}")
    print(f"                 ≈ {float(grad_p_AA):.4f}\n")

    print("Gradient of V w.r.t. pi(A->B|A):")
    print(f"dV/d(pi(A->B|A)) = d_mu(A) * Q_pi(A, A->B)")
    print(f"                 = {f_d_A} * {f_Q_A_to_B}")
    print(f"                 = {f_grad_p_AB}")
    print(f"                 ≈ {float(grad_p_AB):.4f}")

solve_policy_gradient()