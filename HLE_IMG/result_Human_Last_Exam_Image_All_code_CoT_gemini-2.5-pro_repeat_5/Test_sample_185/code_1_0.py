import numpy as np
from fractions import Fraction

def solve_policy_gradient():
    """
    Solves the policy gradient problem as described.
    """
    # Step 1: Define MDP parameters
    gamma = Fraction(1, 2)
    mu0 = np.array([Fraction(1, 4), Fraction(3, 4)])  # [mu0(A), mu0(B)]

    # Policy pi_1
    pi_AA = Fraction(1, 3)
    pi_AB = Fraction(2, 3)
    pi_BA = Fraction(1, 2)
    pi_BB = Fraction(1, 2)

    # Rewards r(s -> s')
    r_AA = Fraction(1)
    r_AB = Fraction(0)
    r_BA = Fraction(1)
    r_BB = Fraction(0)

    # Step 2: Calculate State-Value functions (V_A, V_B)
    # System of equations: (I - gamma * P) * V = R_pi
    # P is the transition matrix under pi_1
    P = np.array([[pi_AA, pi_AB],
                  [pi_BA, pi_BB]])

    # R_pi is the expected reward vector from each state
    R_pi = np.array([pi_AA * r_AA + pi_AB * r_AB,
                     pi_BA * r_BA + pi_BB * r_BB])

    # A * V = R_pi, where A = (I - gamma * P)
    A_V = np.identity(2) - gamma * P
    V = np.linalg.solve(A_V, R_pi)
    V_A, V_B = V[0], V[1]

    print("--- Intermediate Calculations ---")
    print(f"Calculated State-Value V(A) = {V_A}")
    print(f"Calculated State-Value V(B) = {V_B}\n")

    # Step 3: Calculate Action-Value functions (Q-values)
    Q_AA = r_AA + gamma * V_A
    Q_AB = r_AB + gamma * V_B
    # Q_BA and Q_BB are not needed for the gradients wrt pi(A|A)
    # Q_BA = r_BA + gamma * V_A
    # Q_BB = r_BB + gamma * V_B

    print(f"Calculated Action-Value Q(A, A->A) = {Q_AA}")
    print(f"Calculated Action-Value Q(A, A->B) = {Q_AB}\n")

    # Step 4: Calculate discounted state visitation frequencies (d_A, d_B)
    # System of equations: d' * (I - gamma * P) = mu0'  (d' is a row vector)
    # Transposing this gives: (I - gamma * P)' * d = mu0 => (I' - gamma * P') * d = mu0
    # Since I is symmetric: (I - gamma * P.T) * d = mu0
    A_d = np.identity(2) - gamma * P.T
    d = np.linalg.solve(A_d, mu0)
    d_A, d_B = d[0], d[1]
    
    print(f"Calculated Discounted State Visitation Frequency d(A) = {d_A}")
    print(f"Calculated Discounted State Visitation Frequency d(B) = {d_B}\n")

    # Step 5: Apply the Policy Gradient Theorem
    # Since pi(A->B|A) = 1 - pi(A->A|A), we can use the advantage form.
    # dV/d(pi(A->A|A)) = d(A) * (Q(A,A->A) - Q(A,A->B))
    grad_AA = d_A * (Q_AA - Q_AB)
    
    # dV/d(pi(A->B|A)) = d(A) * (Q(A,A->B) - Q(A,A->A))
    grad_AB = d_A * (Q_AB - Q_AA)

    print("--- Final Gradient Calculations ---")
    print("Gradient dV/d(pi(A->A|A)):")
    print(f"d(A) * (Q(A, A->A) - Q(A, A->B)) = {d_A} * ({Q_AA} - {Q_AB})")
    print(f"= {d_A} * ({Q_AA - Q_AB})")
    print(f"= {grad_AA}\n")
    
    print("Gradient dV/d(pi(A->B|A)):")
    print(f"d(A) * (Q(A, A->B) - Q(A, A->A)) = {d_A} * ({Q_AB} - {Q_AA})")
    print(f"= {d_A} * ({Q_AB - Q_AA})")
    print(f"= {grad_AB}")
    
    # Returning the final numerical answers in the required format
    # The problem asks for two gradients, so we provide them both.
    final_answer = f"{grad_AA}, {grad_AB}"
    # This is an unusual format, but we follow it.
    return final_answer

final_answer_string = solve_policy_gradient()
# The required format is <<<answer content>>>. The question asks for two values.
# I will output them as a comma-separated string.
print(f"\n<<<{final_answer_string}>>>")