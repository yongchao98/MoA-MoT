import numpy as np
from fractions import Fraction

def solve_mdp_gradients():
    """
    This function calculates the policy gradients for the given MDP.
    """
    # Step 1: Define MDP and Policy Parameters
    gamma = Fraction(1, 2)
    mu0 = np.array([Fraction(1, 4), Fraction(3, 4)])

    # Policy pi_1
    pi = {
        ('A', 'A'): Fraction(1, 3),
        ('A', 'B'): Fraction(2, 3),
        ('B', 'A'): Fraction(1, 2),
        ('B', 'B'): Fraction(1, 2)
    }

    # Rewards r(s -> s')
    rewards = {
        ('A', 'A'): Fraction(1),
        ('A', 'B'): Fraction(0),
        ('B', 'A'): Fraction(1),
        ('B', 'B'): Fraction(0)
    }

    # Step 2: Calculate State-Value Function V^pi(s)
    # We solve the system of Bellman equations:
    # V(A) = pi(A->A|A)[r(A->A) + gV(A)] + pi(A->B|A)[r(A->B) + gV(B)]
    # V(B) = pi(B->A|B)[r(B->A) + gV(A)] + pi(B->B|B)[r(B->B) + gV(B)]
    #
    # (1 - g*pi(A->A|A))V(A) - g*pi(A->B|A)V(B) = pi(A->A|A)r(A->A) + pi(A->B|A)r(A->B)
    # -g*pi(B->A|B)V(A) + (1 - g*pi(B->B|B))V(B) = pi(B->A|B)r(B->A) + pi(B->B|B)r(B->B)
    
    # Coefficients for the linear system M*v = c
    M = np.array([
        [1 - gamma * pi[('A', 'A')], -gamma * pi[('A', 'B')]],
        [-gamma * pi[('B', 'A')], 1 - gamma * pi[('B', 'B')]]
    ])
    
    c = np.array([
        pi[('A', 'A')] * rewards[('A', 'A')] + pi[('A', 'B')] * rewards[('A', 'B')],
        pi[('B', 'A')] * rewards[('B', 'A')] + pi[('B', 'B')] * rewards[('B', 'B')]
    ])

    # Solve for V = [V(A), V(B)]
    V = np.linalg.solve(M, c)
    V_A, V_B = V[0], V[1]

    print("--- Intermediate Calculations ---")
    print(f"1. State-Value Function:")
    print(f"V(A) = {V_A}")
    print(f"V(B) = {V_B}\n")

    # Step 3: Calculate Action-Value Function Q^pi(s, a)
    # Q(s, s->s') = r(s->s') + gamma * V(s')
    Q_A_A = rewards[('A', 'A')] + gamma * V_A
    Q_A_B = rewards[('A', 'B')] + gamma * V_B
    
    print(f"2. Action-Value Function (for state A):")
    print(f"Q(A, A->A) = r(A->A) + γ*V(A) = {rewards[('A', 'A')]} + {gamma}*{V_A} = {Q_A_A}")
    print(f"Q(A, A->B) = r(A->B) + γ*V(B) = {rewards[('A', 'B')]} + {gamma}*{V_B} = {Q_A_B}\n")

    # Step 4: Calculate Discounted State Visitation Frequency d^pi(s)
    # d^T = mu0^T * (I - gamma * P^pi)^-1
    P_pi = np.array([
        [pi[('A', 'A')], pi[('A', 'B')]],
        [pi[('B', 'A')], pi[('B', 'B')]]
    ])
    
    I = np.identity(2, dtype=object)
    M_d_inv = np.linalg.inv(I - gamma * P_pi)
    
    d = mu0 @ M_d_inv
    d_A, d_B = d[0], d[1]

    print(f"3. Discounted State Visitation Frequency:")
    print(f"d(A) = {d_A}")
    print(f"d(B) = {d_B}\n")

    # Step 5: Compute the Gradients
    # grad(θ_s,a) = d(s) * (Q(s,a) - V(s)) is not quite right for constrained params.
    # The correct derivation for ∂V/∂π(A→A|A) is d(A) * [Q(A, A→A) - Q(A, A→B)]
    grad_pi_A_A = d_A * (Q_A_A - Q_A_B)
    grad_pi_A_B = d_A * (Q_A_B - Q_A_A) # which is -grad_pi_A_A

    print("--- Final Gradient Computations ---")
    print("1. Gradient ∂V/∂π(A→A|A):")
    print("   Formula: d(A) * [Q(A, A→A) - Q(A, A→B)]")
    print(f"   = {d_A} * [{Q_A_A} - {Q_A_B}]")
    print(f"   = {d_A} * [{Q_A_A - Q_A_B}]")
    print(f"   = {grad_pi_A_A}\n")

    print("2. Gradient ∂V/∂π(A→B|A):")
    print("   Formula: d(A) * [Q(A, A→B) - Q(A, A→A)]")
    print(f"   = {d_A} * [{Q_A_B} - {Q_A_A}]")
    print(f"   = {d_A} * [{Q_A_B - Q_A_A}]")
    print(f"   = {grad_pi_A_B}\n")
    
    return grad_pi_A_A, grad_pi_A_B

if __name__ == '__main__':
    grad1, grad2 = solve_mdp_gradients()
    # The final answer format requested by the user
    print(f"<<<{grad1}, {grad2}>>>")
