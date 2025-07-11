import numpy as np
from fractions import Fraction

def solve_policy_gradient():
    """
    Solves for the policy gradients of the given MDP.
    """
    # Step 0: Define the MDP components
    gamma = Fraction(1, 2)
    mu0 = {'A': Fraction(1, 4), 'B': Fraction(3, 4)}
    pi1 = {
        'A': {'A->A': Fraction(1, 3), 'A->B': Fraction(2, 3)},
        'B': {'B->A': Fraction(1, 2), 'B->B': Fraction(1, 2)}
    }
    R = {
        'A': {'A->A': 1, 'A->B': 0},
        'B': {'B->A': 1, 'B->B': 0}
    }
    states = ['A', 'B']

    # --- Step 1: Calculate V^pi(s) ---
    # The Bellman expectation equation gives a system of linear equations:
    # (I - gamma * P_pi) * V = R_pi
    
    # Transition matrix P_pi under policy pi1
    P_pi = np.array([
        [pi1['A']['A->A'], pi1['A']['A->B']],  # P(s'|s=A)
        [pi1['B']['B->A'], pi1['B']['B->B']]   # P(s'|s=B)
    ], dtype=float)

    # Expected immediate reward R_pi from each state
    R_pi_vec = np.array([
        pi1['A']['A->A'] * R['A']['A->A'] + pi1['A']['A->B'] * R['A']['A->B'],
        pi1['B']['B->A'] * R['B']['B->A'] + pi1['B']['B->B'] * R['B']['B->B']
    ], dtype=float)

    # Solve the system (I - gamma * P_pi)V = R_pi for V
    A_mat = np.identity(len(states)) - float(gamma) * P_pi
    V_vec = np.linalg.solve(A_mat, R_pi_vec)
    V = {s: Fraction(v).limit_denominator() for s, v in zip(states, V_vec)}

    # --- Step 2: Calculate Q^pi(s, a) ---
    # Q(s, s->s_next) = R(s, s->s_next) + gamma * V(s_next)
    Q = {
        'A': {
            'A->A': Fraction(R['A']['A->A']) + gamma * V['A'],
            'A->B': Fraction(R['A']['A->B']) + gamma * V['B']
        },
        'B': {
            'B->A': Fraction(R['B']['B->A']) + gamma * V['A'],
            'B->B': Fraction(R['B']['B->B']) + gamma * V['B']
        }
    }

    # --- Step 3: Calculate discounted state visitation frequencies d^pi(s) ---
    # d^T = mu0^T + gamma * d^T * P_pi  =>  d^T * (I - gamma * P_pi) = mu0^T
    # We can solve this by d = mu0 * (I - gamma * P_pi)^-1
    mu0_vec = np.array([float(mu0['A']), float(mu0['B'])])
    d_vec = mu0_vec @ np.linalg.inv(A_mat)
    d = {s: Fraction(dv).limit_denominator() for s, dv in zip(states, d_vec)}

    # --- Step 4: Compute and print the gradients ---
    # Because pi(A->A|A) + pi(A->B|A) = 1, we have d(pi(A->B|A))/d(pi(A->A|A)) = -1.
    # The gradient simplifies to: d(A) * [Q(A,A->A) - Q(A,A->B)]
    grad_AA = d['A'] * (Q['A']['A->A'] - Q['A']['A->B'])
    
    print("The first gradient is ∂V/∂π(A->A|A):")
    print(f"  = d^π(A) * [Q^π(A, A->A) - Q^π(A, A->B)]")
    print(f"  = {d['A']} * [{Q['A']['A->A']} - {Q['A']['A->B']}]")
    print(f"  = {d['A']} * [{Q['A']['A->A'] - Q['A']['A->B']}]")
    print(f"  = {grad_AA}")
    print(f"  ≈ {float(grad_AA):.4f}")
    print("-" * 30)
    
    # The second gradient is the negative of the first.
    grad_AB = d['A'] * (Q['A']['A->B'] - Q['A']['A->A'])

    print("The second gradient is ∂V/∂π(A->B|A):")
    print(f"  = d^π(A) * [Q^π(A, A->B) - Q^π(A, A->A)]")
    print(f"  = {d['A']} * [{Q['A']['A->B']} - {Q['A']['A->A']}]")
    print(f"  = {d['A']} * [{Q['A']['A->B'] - Q['A']['A->A']}]")
    print(f"  = {grad_AB}")
    print(f"  ≈ {float(grad_AB):.4f}")

    return grad_AA, grad_AB

grad1, grad2 = solve_policy_gradient()
# The final answer format is not specified as a tuple, but rather a final output from the program. 
# The prints within the function serve as the main output. 
# Let's provide the tuple of fractional values as the final answer signal.
# <<< (108/169, -108/169) >>>