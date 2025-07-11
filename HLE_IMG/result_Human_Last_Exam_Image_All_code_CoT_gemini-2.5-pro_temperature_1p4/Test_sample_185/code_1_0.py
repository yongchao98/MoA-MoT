import numpy as np
from fractions import Fraction

def solve_mdp_gradients():
    """
    Computes the policy gradients for the given MDP.
    """
    # Step 1: Define MDP parameters from the problem description
    gamma = Fraction(1, 2)
    mu0 = np.array([Fraction(1, 4), Fraction(3, 4)])  # [mu0(A), mu0(B)]

    # Policy pi_1
    pi = {
        'A': {'A': Fraction(1, 3), 'B': Fraction(2, 3)},
        'B': {'A': Fraction(1, 2), 'B': Fraction(1, 2)}
    }

    # Rewards r(s -> s')
    rewards = {
        'A': {'A': 1, 'B': 0},
        'B': {'A': 1, 'B': 0}
    }
    
    print("Step 1: MDP Parameters")
    print(f"gamma = {gamma}")
    print(f"mu0(A) = {mu0[0]}, mu0(B) = {mu0[1]}")
    print(f"pi(A->A|A) = {pi['A']['A']}, pi(A->B|A) = {pi['A']['B']}")
    print(f"pi(B->A|B) = {pi['B']['A']}, pi(B->B|B) = {pi['B']['B']}")
    print("-" * 30)

    # Step 2: Calculate State Values V(s)
    # The Bellman equations are V = R_pi + gamma * P_pi * V
    # which can be rearranged to (I - gamma * P_pi) * V = R_pi
    
    # State transition matrix P_pi
    P_pi = np.array([
        [pi['A']['A'], pi['A']['B']],
        [pi['B']['A'], pi['B']['B']]
    ])

    # Expected immediate rewards vector R_pi
    R_pi = np.array([
        pi['A']['A'] * rewards['A']['A'] + pi['A']['B'] * rewards['A']['B'],
        pi['B']['A'] * rewards['B']['A'] + pi['B']['B'] * rewards['B']['B']
    ])
    
    I = np.identity(2, dtype=Fraction)
    A_mat = I - gamma * P_pi
    
    # Solve for V
    V_pi = np.linalg.solve(A_mat, R_pi)
    V_A, V_B = V_pi[0], V_pi[1]

    print("Step 2: State Values V(s)")
    print(f"V(A) = {V_A}")
    print(f"V(B) = {V_B}")
    print("-" * 30)
    
    # Step 3: Calculate Action-Values Q(s, a)
    # Q(s, s') = r(s->s') + gamma * V(s')
    Q_AA = rewards['A']['A'] + gamma * V_A
    Q_AB = rewards['A']['B'] + gamma * V_B
    Q_BA = rewards['B']['A'] + gamma * V_A
    Q_BB = rewards['B']['B'] + gamma * V_B

    print("Step 3: Action-Values Q(s, a)")
    print(f"Q(A, A->A) = {Q_AA}")
    print(f"Q(A, A->B) = {Q_AB}")
    print(f"Q(B, B->A) = {Q_BA}")
    print(f"Q(B, B->B) = {Q_BB}")
    print("-" * 30)

    # Step 4: Calculate Discounted State Visitation Frequencies d(s)
    # d_T = mu0_T * (I - gamma * P_pi)^-1
    inv_A_mat = np.linalg.inv(A_mat)
    d_pi = mu0.dot(inv_A_mat)
    d_A, d_B = d_pi[0], d_pi[1]

    print("Step 4: Discounted State Visitation Frequencies d(s)")
    print(f"d(A) = {d_A}")
    print(f"d(B) = {d_B}")
    print("-" * 30)

    # Step 5: Apply the Policy Gradient Theorem
    # For a policy parameterized by theta_AA = pi(A->A|A), with pi(A->B|A) = 1 - theta_AA,
    # the gradient is: d(A) * [ (1)*Q(A, A->A) + (-1)*Q(A, A->B) ]
    grad_V_theta_AA = d_A * (Q_AA - Q_AB)

    # For a policy parameterized by theta_AB = pi(A->B|A), with pi(A->A|A) = 1 - theta_AB,
    # the gradient is: d(A) * [ (-1)*Q(A, A->A) + (1)*Q(A, A->B) ]
    grad_V_theta_AB = d_A * (Q_AB - Q_AA)
    
    print("Step 5: Compute Gradients\n")
    
    print("Gradient of V w.r.t. pi(A->A|A):")
    print(f"Formula: d(A) * [Q(A, A->A) - Q(A, A->B)]")
    print(f"= {d_A} * [{Q_AA} - {Q_AB}]")
    print(f"= {d_A} * [{Q_AA - Q_AB}]")
    print(f"= {grad_V_theta_AA}\n")

    print("Gradient of V w.r.t. pi(A->B|A):")
    print(f"Formula: d(A) * [Q(A, A->B) - Q(A, A->A)]")
    print(f"= {d_A} * [{Q_AB} - {Q_AA}]")
    print(f"= {d_A} * [{Q_AB - Q_AA}]")
    print(f"= {grad_V_theta_AB}\n")
    
    # Final answer
    return grad_V_theta_AA, grad_V_theta_AB

if __name__ == '__main__':
    grad1, grad2 = solve_mdp_gradients()
    print(f"<<<{grad1}, {grad2}>>>")
