from fractions import Fraction

def solve_mdp_gradients():
    """
    Solves for the policy gradients of a given 2-state MDP.
    """
    # 1. Define MDP parameters
    gamma = Fraction(1, 2)
    mu0 = [Fraction(1, 4), Fraction(3, 4)]  # [mu(A), mu(B)]
    states = ['A', 'B']
    s_map = {s: i for i, s in enumerate(states)}

    # Policy pi_1
    pi = [
        [Fraction(1, 3), Fraction(2, 3)],  # From A: P(A->A), P(A->B)
        [Fraction(1, 2), Fraction(1, 2)]   # From B: P(B->A), P(B->B)
    ]

    # Rewards r(s, s')
    R = [
        [Fraction(1), Fraction(0)],  # r(A,A), r(A,B)
        [Fraction(1), Fraction(0)]   # r(B,A), r(B,B)
    ]

    # 2. Calculate V^pi
    # V^pi = (I - gamma * P_pi)^-1 * R_pi
    # where R_pi(s) = E[r | s] = sum_s' pi(s'|s) * r(s, s')
    P_pi = pi
    R_pi = [sum(pi[i][j] * R[i][j] for j in range(len(states))) for i in range(len(states))]

    # Build the matrix M = (I - gamma * P_pi)
    M = [
        [Fraction(1) - gamma * P_pi[0][0], -gamma * P_pi[0][1]],
        [-gamma * P_pi[1][0], Fraction(1) - gamma * P_pi[1][1]]
    ]

    # Invert the 2x2 matrix M
    det_M = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    M_inv = [
        [M[1][1] / det_M, -M[0][1] / det_M],
        [-M[1][0] / det_M, M[0][0] / det_M]
    ]

    # Solve for V_pi = M_inv * R_pi
    V_pi = [
        M_inv[0][0] * R_pi[0] + M_inv[0][1] * R_pi[1],
        M_inv[1][0] * R_pi[0] + M_inv[1][1] * R_pi[1]
    ]
    V_A, V_B = V_pi[s_map['A']], V_pi[s_map['B']]

    # 3. Calculate Q^pi
    # Q(s, s') = r(s, s') + gamma * V(s')
    Q_A_to_A = R[s_map['A']][s_map['A']] + gamma * V_A
    Q_A_to_B = R[s_map['A']][s_map['B']] + gamma * V_B

    # 4. Calculate discounted state visitation frequencies d^pi
    # d_pi = mu0^T * M_inv
    d_pi = [
        mu0[0] * M_inv[0][0] + mu0[1] * M_inv[1][0],
        mu0[0] * M_inv[0][1] + mu0[1] * M_inv[1][1]
    ]
    d_A = d_pi[s_map['A']]

    # 5. Calculate gradients
    grad_A_to_A = d_A * (Q_A_to_A - Q_A_to_B)
    grad_A_to_B = d_A * (Q_A_to_B - Q_A_to_A)

    # 6. Print the results step-by-step
    print("--- Intermediate Values ---")
    print(f"V(A) = {V_A}, V(B) = {V_B}")
    print(f"Q(A, A->A) = {Q_A_to_A}")
    print(f"Q(A, A->B) = {Q_A_to_B}")
    print(f"d(A) = {d_A}")
    print("\n--- Gradient Calculations ---")

    print("\nGradient of V w.r.t. π(A->A|A):")
    print("∂V^π(μ₀)/∂π(A->A|A) = d^π(A) * (Q^π(A, A->A) - Q^π(A, A->B))")
    print(f"= {d_A} * ({Q_A_to_A} - {Q_A_to_B})")
    print(f"= {d_A} * ({Q_A_to_A - Q_A_to_B})")
    print(f"= {grad_A_to_A}")

    print("\nGradient of V w.r.t. π(A->B|A):")
    print("∂V^π(μ₀)/∂π(A->B|A) = d^π(A) * (Q^π(A, A->B) - Q^π(A, A->A))")
    print(f"= {d_A} * ({Q_A_to_B} - {Q_A_to_A})")
    print(f"= {d_A} * ({Q_A_to_B - Q_A_to_A})")
    print(f"= {grad_A_to_B}")

if __name__ == '__main__':
    solve_mdp_gradients()