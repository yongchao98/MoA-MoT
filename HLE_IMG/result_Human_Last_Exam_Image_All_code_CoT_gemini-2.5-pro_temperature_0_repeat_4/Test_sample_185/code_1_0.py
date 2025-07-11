import numpy as np
from fractions import Fraction

def solve_policy_gradient():
    """
    Solves for the policy gradients for the given MDP.
    """
    # --- 1. Define MDP parameters ---
    gamma = 0.5
    mu0 = {'A': 0.25, 'B': 0.75}
    pi = {
        ('A', 'A'): 1/3,
        ('A', 'B'): 2/3,
        ('B', 'A'): 1/2,
        ('B', 'B'): 1/2
    }
    r = {
        ('A', 'A'): 1,
        ('A', 'B'): 0,
        ('B', 'A'): 1,
        ('B', 'B'): 0
    }
    states = ['A', 'B']

    # --- 2. Calculate State-Value Function V(s) ---
    # The Bellman equations form a linear system M_v * V = C_v
    # (1 - g*pi(A->A|A))V_A - g*pi(A->B|A)V_B = pi(A->A|A)r(A->A) + pi(A->B|A)r(A->B)
    # -g*pi(B->A|B)V_A + (1 - g*pi(B->B|B))V_B = pi(B->A|B)r(B->A) + pi(B->B|B)r(B->B)
    M_v = np.array([
        [1 - gamma * pi[('A', 'A')], -gamma * pi[('A', 'B')]],
        [-gamma * pi[('B', 'A')], 1 - gamma * pi[('B', 'B')]]
    ])
    C_v = np.array([
        pi[('A', 'A')] * r[('A', 'A')] + pi[('A', 'B')] * r[('A', 'B')],
        pi[('B', 'A')] * r[('B', 'A')] + pi[('B', 'B')] * r[('B', 'B')]
    ])
    V = np.linalg.solve(M_v, C_v)
    V_A, V_B = V[0], V[1]

    # --- 3. Calculate Action-Value Function Q(s, a) ---
    # Q(s, a) = r(s, a) + gamma * V(s')
    Q_A_to_A = r[('A', 'A')] + gamma * V_A
    Q_A_to_B = r[('A', 'B')] + gamma * V_B

    # --- 4. Calculate Discounted State Visitation Frequency d(s) ---
    # The equations for d(s) form a linear system M_d * D = C_d
    # (1 - g*pi(A->A|A))d_A - g*pi(B->A|B)d_B = mu0(A)
    # -g*pi(A->B|A)d_A + (1 - g*pi(B->B|B))d_B = mu0(B)
    M_d = np.array([
        [1 - gamma * pi[('A', 'A')], -gamma * pi[('B', 'A')]],
        [-gamma * pi[('A', 'B')], 1 - gamma * pi[('B', 'B')]]
    ])
    C_d = np.array([mu0['A'], mu0['B']])
    D = np.linalg.solve(M_d, C_d)
    d_A, d_B = D[0], D[1]

    # --- 5. Compute the Gradients ---
    grad_pi_A_to_A = d_A * Q_A_to_A
    grad_pi_A_to_B = d_A * Q_A_to_B

    # --- 6. Print results using fractions for exactness ---
    V_A_f = Fraction(V_A).limit_denominator()
    V_B_f = Fraction(V_B).limit_denominator()
    Q_A_to_A_f = Fraction(Q_A_to_A).limit_denominator()
    Q_A_to_B_f = Fraction(Q_A_to_B).limit_denominator()
    d_A_f = Fraction(d_A).limit_denominator()
    grad1_f = Fraction(grad_pi_A_to_A).limit_denominator()
    grad2_f = Fraction(grad_pi_A_to_B).limit_denominator()

    print("--- Intermediate Calculations ---")
    print(f"State-Value V(A) = {V_A_f.numerator}/{V_A_f.denominator}")
    print(f"State-Value V(B) = {V_B_f.numerator}/{V_B_f.denominator}")
    print(f"Discounted State Visitation d(A) = {d_A_f.numerator}/{d_A_f.denominator}")
    print("-" * 35)
    print("\n--- Final Gradient Computations ---")

    # Gradient 1
    print(f"1. Gradient for π(A->A|A):")
    print(f"   ∂V/∂π(A->A|A) = d(A) * Q(A, A->A)")
    print(f"   Q(A, A->A) = r(A->A) + γ*V(A) = {r[('A', 'A')]} + {gamma}*({V_A_f.numerator}/{V_A_f.denominator}) = {Q_A_to_A_f.numerator}/{Q_A_to_A_f.denominator}")
    print(f"   ∂V/∂π(A->A|A) = ({d_A_f.numerator}/{d_A_f.denominator}) * ({Q_A_to_A_f.numerator}/{Q_A_to_A_f.denominator}) = {grad1_f.numerator}/{grad1_f.denominator}")
    print(f"   Result ≈ {grad_pi_A_to_A:.4f}\n")

    # Gradient 2
    print(f"2. Gradient for π(A->B|A):")
    print(f"   ∂V/∂π(A->B|A) = d(A) * Q(A, A->B)")
    print(f"   Q(A, A->B) = r(A->B) + γ*V(B) = {r[('A', 'B')]} + {gamma}*({V_B_f.numerator}/{V_B_f.denominator}) = {Q_A_to_B_f.numerator}/{Q_A_to_B_f.denominator}")
    print(f"   ∂V/∂π(A->B|A) = ({d_A_f.numerator}/{d_A_f.denominator}) * ({Q_A_to_B_f.numerator}/{Q_A_to_B_f.denominator}) = {grad2_f.numerator}/{grad2_f.denominator}")
    print(f"   Result ≈ {grad_pi_A_to_B:.4f}")

if __name__ == '__main__':
    solve_policy_gradient()