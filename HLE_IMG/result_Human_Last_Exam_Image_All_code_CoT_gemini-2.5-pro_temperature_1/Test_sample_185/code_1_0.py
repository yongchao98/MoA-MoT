import numpy as np
from fractions import Fraction

def solve_policy_gradient():
    """
    Solves for the policy gradients of a given MDP.
    """
    # Define the parameters of the MDP and policy
    gamma = Fraction(1, 2)

    # Policy pi_1
    p_AA = Fraction(1, 3)
    p_AB = Fraction(2, 3)
    p_BA = Fraction(1, 2)
    p_BB = Fraction(1, 2)

    # Rewards from the figure
    R_AA = 1
    R_AB = 0
    R_BA = 1
    R_BB = 0

    # Initial distribution
    mu_A = Fraction(1, 4)
    mu_B = Fraction(3, 4)

    # --- Step 1: Calculate the value functions V(A) and V(B) ---
    # The Bellman equations form a linear system A*V = b:
    # (1 - p_AA*gamma) V_A - p_AB*gamma V_B = p_AA*R_AA + p_AB*R_AB
    # -p_BA*gamma V_A + (1 - p_BB*gamma) V_B = p_BA*R_BA + p_BB*R_BB
    A_matrix = np.array([
        [float(1 - p_AA * gamma), float(-p_AB * gamma)],
        [float(-p_BA * gamma), float(1 - p_BB * gamma)]
    ])
    b_vector = np.array([
        float(p_AA * R_AA + p_AB * R_AB),
        float(p_BA * R_BA + p_BB * R_BB)
    ])
    V = np.linalg.solve(A_matrix, b_vector)
    V_A, V_B = Fraction(V[0]).limit_denominator(), Fraction(V[1]).limit_denominator()

    # --- Step 2: Calculate the derivatives of the value functions w.r.t p = p_AA ---
    # Differentiating the Bellman equations w.r.t. p gives a linear system C*V' = d
    p = p_AA
    C_matrix = np.array([
        [float(1 - p * gamma), float(-(1 - p) * gamma)],
        [float(-p_BA * gamma), float(1 - p_BB * gamma)]
    ])
    d_vector = np.array([
        float(R_AA + gamma * V_A - gamma * V_B),
        0
    ])
    V_prime = np.linalg.solve(C_matrix, d_vector)
    V_prime_A, V_prime_B = Fraction(V_prime[0]).limit_denominator(), Fraction(V_prime[1]).limit_denominator()

    # --- Step 3: Compute the final gradients ---
    grad_p_AA = mu_A * V_prime_A + mu_B * V_prime_B
    grad_p_AB = -grad_p_AA

    # --- Step 4: Print the results and calculations ---
    print("Let p = pi(A->A|A) = 1/3. Then pi(A->B|A) = 1 - p = 2/3.")
    print(f"First, we compute the state values V(A) and V(B), which are {V_A} and {V_B} respectively.")
    print(f"Next, by differentiating the Bellman equations, we find the derivatives dV(A)/dp = {V_prime_A} and dV(B)/dp = {V_prime_B}.")
    
    print("\n1. To compute the gradient dV/d(pi(A->A|A)):")
    term1 = mu_A * V_prime_A
    term2 = mu_B * V_prime_B
    print(f"   dV/d(pi(A->A|A)) = mu_0(A) * dV(A)/dp + mu_0(B) * dV(B)/dp")
    print(f"                   = ({mu_A}) * ({V_prime_A}) + ({mu_B}) * ({V_prime_B})")
    print(f"                   = {term1} + {term2}")
    print(f"                   = {grad_p_AA}")

    print("\n2. To compute the gradient dV/d(pi(A->B|A)):")
    print(f"   dV/d(pi(A->B|A)) = -dV/d(pi(A->A|A))")
    print(f"                   = {grad_p_AB}")

solve_policy_gradient()