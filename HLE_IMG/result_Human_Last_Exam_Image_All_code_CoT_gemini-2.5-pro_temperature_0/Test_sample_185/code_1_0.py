import numpy as np
from fractions import Fraction

def solve_mdp_gradients():
    """
    Solves for the policy gradients of a given MDP.
    """
    # --- 1. Define MDP and Policy Parameters ---
    # Policy pi_1
    pi = {
        ('A', 'A'): Fraction(1, 3),
        ('A', 'B'): Fraction(2, 3),
        ('B', 'A'): Fraction(1, 2),
        ('B', 'B'): Fraction(1, 2)
    }

    # Rewards r(s -> s')
    R = {
        ('A', 'A'): 1,
        ('A', 'B'): 0,
        ('B', 'A'): 1,
        ('B', 'B'): 0
    }

    # Discount factor
    gamma = Fraction(1, 2)

    # Initial distribution mu_0
    mu0 = {'A': Fraction(1, 4), 'B': Fraction(3, 4)}

    states = ['A', 'B']

    # --- 2. Calculate State-Value Functions V(s) ---
    # We solve the Bellman equation system: (I - gamma * P_pi) * V = R_pi
    print("Step 1: Calculate State-Value Functions V(s)")
    
    # Transition matrix P_pi
    P_pi = np.array([
        [pi[('A', 'A')], pi[('A', 'B')]],  # Transitions from A
        [pi[('B', 'A')], pi[('B', 'B')]]   # Transitions from B
    ], dtype=object)

    # Expected immediate reward vector R_pi
    R_pi = np.array([
        pi[('A', 'A')] * R[('A', 'A')] + pi[('A', 'B')] * R[('A', 'B')],
        pi[('B', 'A')] * R[('B', 'A')] + pi[('B', 'B')] * R[('B', 'B')]
    ], dtype=object)

    # Identity matrix
    I = np.identity(len(states), dtype=object)

    # Coefficient matrix for the linear system
    A_matrix = I - gamma * P_pi

    # Solve for V
    V_vec = np.linalg.solve(A_matrix.astype(float), R_pi.astype(float))
    V = {s: Fraction(v).limit_denominator() for s, v in zip(states, V_vec)}

    print(f"Solving the Bellman equations for policy pi_1 yields:")
    print(f"V(A) = {V['A']}")
    print(f"V(B) = {V['B']}")
    print("-" * 40)

    # --- 3. Calculate Action-Value Functions Q(s, a) ---
    print("Step 2: Calculate Action-Value Functions Q(s, a)")
    Q = {}
    Q[('A', 'A')] = R[('A', 'A')] + gamma * V['A']
    Q[('A', 'B')] = R[('A', 'B')] + gamma * V['B']
    
    print(f"Using Q(s, s->s') = r(s->s') + gamma * V(s'):")
    print(f"Q(A, A->A) = {R[('A', 'A')]} + {gamma} * {V['A']} = {Q[('A', 'A')]}")
    print(f"Q(A, A->B) = {R[('A', 'B')]} + {gamma} * {V['B']} = {Q[('A', 'B')]}")
    print("-" * 40)

    # --- 4. Calculate Discounted State Visitation Frequencies d(s) ---
    print("Step 3: Calculate Discounted State Visitation Frequency d(A)")
    mu0_vec = np.array([mu0['A'], mu0['B']], dtype=object)
    A_inv = np.linalg.inv(A_matrix.astype(float))
    d_vec = mu0_vec.astype(float) @ A_inv
    d = {s: Fraction(val).limit_denominator() for s, val in zip(states, d_vec)}

    print("Using d^T = mu0^T * (I - gamma * P_pi)^-1:")
    print(f"d(A) = {d['A']}")
    print("-" * 40)

    # --- 5. Compute the Gradients ---
    print("Step 4: Compute the Gradients")
    grad_A_to_A = d['A'] * (Q[('A', 'A')] - Q[('A', 'B')])
    grad_A_to_B = -grad_A_to_A

    print("The first gradient is calculated as: dV/d(pi(A->A|A)) = d(A) * [Q(A, A->A) - Q(A, A->B)]")
    print(f"Result = {d['A']} * [{Q[('A', 'A')]} - {Q[('A', 'B')]}]")
    print(f"       = {d['A']} * [{Q[('A', 'A')] - Q[('A', 'B')]}]")
    print(f"       = {grad_A_to_A}")
    print("\nThe second gradient is the negative of the first due to the policy constraint:")
    print(f"dV/d(pi(A->B|A)) = {grad_A_to_B}")
    
    return grad_A_to_A, grad_A_to_B

if __name__ == '__main__':
    grad1, grad2 = solve_mdp_gradients()
    print(f"\n<<<The two gradients are {grad1} and {grad2}>>>")
