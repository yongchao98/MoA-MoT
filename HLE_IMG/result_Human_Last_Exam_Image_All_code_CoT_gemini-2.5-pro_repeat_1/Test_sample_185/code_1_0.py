import numpy as np
from fractions import Fraction

def solve_policy_gradient():
    """
    Solves the policy gradient problem for the given MDP.
    """
    # Step 1: Define MDP parameters from the problem description
    gamma = Fraction(1, 2)
    mu0 = np.array([Fraction(1, 4), Fraction(3, 4)])

    # Policy pi_1 probabilities
    # pi[state_idx, action_idx] where A=0, B=1
    # For state A, actions are A->A (idx 0), A->B (idx 1)
    # For state B, actions are B->A (idx 0), B->B (idx 1)
    pi_probs = {
        'A': {'A->A': Fraction(1, 3), 'A->B': Fraction(2, 3)},
        'B': {'B->A': Fraction(1, 2), 'B->B': Fraction(1, 2)}
    }

    # Rewards r(s, action)
    rewards = {
        'A': {'A->A': Fraction(1), 'A->B': Fraction(0)},
        'B': {'B->A': Fraction(1), 'B->B': Fraction(0)}
    }

    # Transition dynamics P(s' | s, action) are deterministic
    # P_matrix[s_idx, action_idx, s'_idx]
    # For s=A, a=A->A, s'=A is 1, s'=B is 0.
    # For s=A, a=A->B, s'=A is 0, s'=B is 1.
    # etc.

    # Step 2: Calculate State Values (V_pi)
    # First, we need the transition matrix P_pi and reward vector R_pi under the policy
    # P_pi[i, j] = sum_a pi(a|i) * P(j|i,a)
    P_pi = np.array([
        [pi_probs['A']['A->A'], pi_probs['A']['A->B']],  # From state A to A, B
        [pi_probs['B']['B->A'], pi_probs['B']['B->B']]   # From state B to A, B
    ], dtype=object)

    # R_pi[i] = sum_a pi(a|i) * r(i,a)
    R_pi = np.array([
        pi_probs['A']['A->A'] * rewards['A']['A->A'] + pi_probs['A']['A->B'] * rewards['A']['A->B'],
        pi_probs['B']['B->A'] * rewards['B']['B->A'] + pi_probs['B']['B->B'] * rewards['B']['B->B']
    ], dtype=object)

    # Solve the Bellman equation: (I - gamma * P_pi) * V = R_pi
    I = np.identity(2, dtype=object)
    A_matrix = I - gamma * P_pi
    
    # Invert the 2x2 matrix A_matrix to solve for V
    a, b = A_matrix[0, 0], A_matrix[0, 1]
    c, d = A_matrix[1, 0], A_matrix[1, 1]
    det = a * d - b * c
    A_inv = (1/det) * np.array([[d, -b], [-c, a]], dtype=object)

    V_pi = A_inv @ R_pi
    V_A, V_B = V_pi[0], V_pi[1]

    # Step 3: Calculate Action Values (Q-values) for actions from state A
    # Q(s, action) = r(s, action) + gamma * V(s_next)
    Q_A_AtoA = rewards['A']['A->A'] + gamma * V_A
    Q_A_AtoB = rewards['A']['A->B'] + gamma * V_B

    # Step 4: Calculate Discounted State Visitation Frequencies (d_pi)
    # d_pi^T = mu0^T @ (I - gamma * P_pi)^-1
    d_pi = mu0 @ A_inv
    d_A, d_B = d_pi[0], d_pi[1]

    # Step 5: Compute the gradients
    grad_A_AtoA = d_A * Q_A_AtoA
    grad_A_AtoB = d_A * Q_A_AtoB

    # Print the results and intermediate steps
    print("Policy Gradient Calculation:")
    print("-" * 40)
    print(f"State Values:\n V(A) = {V_A.numerator}/{V_A.denominator}\n V(B) = {V_B.numerator}/{V_B.denominator}")
    print("-" * 40)
    print(f"Discounted State Visitations:\n d(A) = {d_A.numerator}/{d_A.denominator}\n d(B) = {d_B.numerator}/{d_B.denominator}")
    print("-" * 40)
    print(f"Action Values from State A:\n Q(A, A->A) = {Q_A_AtoA.numerator}/{Q_A_AtoA.denominator}\n Q(A, A->B) = {Q_A_AtoB.numerator}/{Q_A_AtoB.denominator}")
    print("-" * 40)

    # Final Gradient 1
    print("Gradient for pi(A->A|A):")
    print(f"∂V/∂π(A->A|A) = d(A) * Q(A, A->A)")
    print(f"             = ({d_A.numerator}/{d_A.denominator}) * ({Q_A_AtoA.numerator}/{Q_A_AtoA.denominator})")
    print(f"             = {grad_A_AtoA.numerator}/{grad_A_AtoA.denominator}")
    print()

    # Final Gradient 2
    print("Gradient for pi(A->B|A):")
    print(f"∂V/∂π(A->B|A) = d(A) * Q(A, A->B)")
    print(f"             = ({d_A.numerator}/{d_A.denominator}) * ({Q_A_AtoB.numerator}/{Q_A_AtoB.denominator})")
    print(f"             = {grad_A_AtoB.numerator}/{grad_A_AtoB.denominator}")

if __name__ == '__main__':
    solve_policy_gradient()