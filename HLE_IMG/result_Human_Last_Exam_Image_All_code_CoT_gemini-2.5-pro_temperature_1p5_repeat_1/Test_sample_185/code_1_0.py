import numpy as np

def calculate_policy_gradients():
    """
    Calculates the policy gradients for the given MDP and policy.
    """
    # Step 0: Define MDP and Policy Parameters
    gamma = 0.5
    mu0 = np.array([0.25, 0.75])  # [mu0(A), mu0(B)]

    # Policy pi_1 probabilities
    pi_AA = 1/3
    pi_AB = 2/3
    pi_BA = 1/2
    pi_BB = 1/2

    # Transition matrix under policy pi_1
    P_pi = np.array([[pi_AA, pi_AB],
                     [pi_BA, pi_BB]])

    # Expected one-step reward vector from each state under pi_1
    R_pi = np.array([
        pi_AA * 1 + pi_AB * 0,  # Expected reward from state A
        pi_BA * 1 + pi_BB * 0   # Expected reward from state B
    ])

    print("Step-by-step calculation of the policy gradients:")

    # Step 1: Calculate State Value Function V^pi
    # We solve the Bellman system: (I - gamma * P_pi) V = R_pi
    I = np.identity(2)
    A_V_system = I - gamma * P_pi
    V = np.linalg.solve(A_V_system, R_pi)
    V_A, V_B = V[0], V[1]
    
    print("\n--- 1. Solving for State Values V^pi(s) ---")
    print("Equation: (I - gamma * P_pi) * V = R_pi")
    print(f"gamma={gamma}, P_pi=\n{P_pi}")
    print(f"Matrix (I - gamma * P_pi):\n{A_V_system}")
    print(f"Vector R_pi (expected rewards): {R_pi}")
    # Manually calculated exact fractions for clarity
    v_a_frac = "10/13"
    v_b_frac = "12/13"
    print(f"Result: V^pi(A) = {V_A:.4f} (exact: {v_a_frac}), V^pi(B) = {V_B:.4f} (exact: {v_b_frac})")


    # Step 2: Calculate Action Value Function Q^pi
    Q_A_A = 1 + gamma * V_A  # Q(A, A->A) = r(A->A) + gamma*V(A)
    Q_A_B = 0 + gamma * V_B  # Q(A, A->B) = r(A->B) + gamma*V(B)
    
    print("\n--- 2. Calculating Action Values Q^pi(A, a) ---")
    # Manually calculated exact fractions
    q_aa_frac = "18/13"
    q_ab_frac = "6/13"
    print(f"Q^pi(A, A->A) = r(A->A) + gamma * V^pi(A) = 1 + {gamma} * {V_A:.4f} = {Q_A_A:.4f} (exact: {q_aa_frac})")
    print(f"Q^pi(A, A->B) = r(A->B) + gamma * V^pi(B) = 0 + {gamma} * {V_B:.4f} = {Q_A_B:.4f} (exact: {q_ab_frac})")

    # Step 3: Calculate Discounted State Visitation Frequency d^pi
    # We solve d * (I - gamma*P_pi) = mu0, which is equivalent to (I - gamma*P_pi)^T * d^T = mu0^T
    A_d_system = A_V_system.T
    d = np.linalg.solve(A_d_system, mu0)
    d_A, d_B = d[0], d[1]
    
    print("\n--- 3. Solving for Discounted State Visitations d^pi(s) ---")
    print("Equation: d * (I - gamma * P_pi) = mu0")
    print(f"mu0 = {mu0}")
    # Manually calculated exact fractions
    d_a_frac = "9/13"
    d_b_frac = "17/13"
    print(f"Result: d^pi(A) = {d_A:.4f} (exact: {d_a_frac}), d^pi(B) = {d_B:.4f} (exact: {d_b_frac})")

    # Step 4: Compute Gradients
    grad_pAA = d_A * (Q_A_A - Q_A_B)
    grad_pAB = d_A * (Q_A_B - Q_A_A)

    print("\n--- 4. Computing the Final Gradients ---")
    print("\nGradient of V w.r.t. pi(A->A|A):")
    print("Formula: d^pi(A) * (Q^pi(A, A->A) - Q^pi(A, A->B))")
    print(f"Calculation: ({d_a_frac}) * (({q_aa_frac}) - ({q_ab_frac}))")
    print(f"             = ({d_a_frac}) * (12/13) = 108/169")
    print(f"Final Value: {grad_pAA:.5f}")

    print("\nGradient of V w.r.t. pi(A->B|A):")
    print("Formula: d^pi(A) * (Q^pi(A, A->B) - Q^pi(A, A->A))")
    print(f"Calculation: ({d_a_frac}) * (({q_ab_frac}) - ({q_aa_frac}))")
    print(f"             = ({d_a_frac}) * (-12/13) = -108/169")
    print(f"Final Value: {grad_pAB:.5f}")

if __name__ == '__main__':
    calculate_policy_gradients()