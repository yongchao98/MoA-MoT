import numpy as np
from fractions import Fraction

def solve_mdp_gradients():
    """
    This function calculates the policy gradients for the given MDP.
    """
    # 1. Define MDP parameters from the problem description
    gamma = Fraction(1, 2)
    # Initial state distribution mu_0(A), mu_0(B)
    mu0 = np.array([Fraction(1, 4), Fraction(3, 4)])

    # Policy pi_1 probabilities
    pi_AA = Fraction(1, 3)
    pi_AB = Fraction(2, 3)
    pi_BA = Fraction(1, 2)
    pi_BB = Fraction(1, 2)

    # Reward function r(s -> s')
    r_AA = 1
    r_AB = 0
    r_BA = 1
    r_BB = 0

    # 2. Calculate State-Value functions V(s) by solving the Bellman equations
    # The Bellman equations in matrix form are: (I - gamma * P) * V = R_pi
    # where P is the transition matrix under the policy and R_pi is the expected reward vector.

    # Transition matrix P
    P = np.array([[pi_AA, pi_AB],
                  [pi_BA, pi_BB]])

    # Expected reward vector R_pi
    R_pi = np.array([pi_AA * r_AA + pi_AB * r_AB,
                     pi_BA * r_BA + pi_BB * r_BB])

    I = np.identity(2)
    # Coefficient matrix for the linear system to find V
    A_V = I - gamma * P

    # Solve for V. numpy.linalg.solve requires floating point numbers.
    V_float = np.linalg.solve(A_V.astype(float), R_pi.astype(float))
    # Convert results back to fractions for precision
    V_A = Fraction(V_float[0]).limit_denominator()
    V_B = Fraction(V_float[1]).limit_denominator()

    print("--- Intermediate Calculations ---")
    print("\nStep a: State-Value Functions V(s)")
    print(f"By solving the Bellman equations, we get:")
    print(f"V(A) = {V_A}")
    print(f"V(B) = {V_B}")

    # 3. Calculate the required Action-Value functions Q(s, a)
    # Q(s, s->s') = r(s->s') + gamma * V(s')
    Q_A_to_A = r_AA + gamma * V_A
    Q_A_to_B = r_AB + gamma * V_B
    
    print("\nStep b: Action-Value Functions Q(s, a)")
    print(f"Q(A, A->A) = r(A->A) + gamma * V(A) = {r_AA} + {gamma} * {V_A} = {Q_A_to_A}")
    print(f"Q(A, A->B) = r(A->B) + gamma * V(B) = {r_AB} + {gamma} * {V_B} = {Q_A_to_B}")

    # 4. Calculate discounted state visitation frequencies d(s)
    # The system is d^T * (I - gamma * P) = mu0^T, which is equivalent to (I - gamma * P)^T * d = mu0.
    A_d = (I - gamma * P).T
    b_d = mu0

    # Solve for d
    d_float = np.linalg.solve(A_d.astype(float), b_d.astype(float))
    d_A = Fraction(d_float[0]).limit_denominator()
    d_B = Fraction(d_float[1]).limit_denominator()

    print("\nStep c: Discounted State Visitation Frequencies d(s)")
    print(f"d(A) = {d_A}")
    print(f"d(B) = {d_B}")

    # 5. Calculate the final gradients
    grad_AA = d_A * Q_A_to_A
    grad_AB = d_A * Q_A_to_B
    
    print("\n--- Final Gradient Computations ---")
    print("\nGradient 1: dV/d(pi(A->A|A)) = d(A) * Q(A, A->A)")
    print(f"= {d_A} * {Q_A_to_A} = {grad_AA}")
    print(f"As a decimal: {float(grad_AA):.4f}")


    print("\nGradient 2: dV/d(pi(A->B|A)) = d(A) * Q(A, A->B)")
    print(f"= {d_A} * {Q_A_to_B} = {grad_AB}")
    print(f"As a decimal: {float(grad_AB):.4f}")

    # Return the final computed values for the answer tag
    return grad_AA, grad_AB

if __name__ == '__main__':
    grad1, grad2 = solve_mdp_gradients()
    # The final answer will be extracted from the print output, but for clarity:
    # print(f"\nFinal Answer: {grad1}, {grad2}")