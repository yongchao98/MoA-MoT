import numpy as np
from fractions import Fraction

def solve_policy_gradient():
    """
    Solves the policy gradient problem for the given MDP.
    """
    # 1. Define MDP and Policy Parameters
    # Policy pi_1
    pi_A_A = 1/3
    pi_A_B = 2/3
    pi_B_A = 1/2
    pi_B_B = 1/2

    # Rewards r(s -> s')
    r_A_A = 1
    r_A_B = 0
    r_B_A = 1
    r_B_B = 0

    # Discount factor
    gamma = 1/2

    # Initial distribution mu_0(s)
    mu0 = np.array([1/4, 3/4]) # [mu(A), mu(B)]

    # 2. Calculate State Values V^pi
    # The Bellman equation for a policy pi is V = R_pi + gamma * P_pi * V,
    # which can be rewritten as (I - gamma * P_pi) * V = R_pi.
    
    # State transition matrix P_pi[i, j] = P(s'=j | s=i, pi)
    P_pi = np.array([
        [pi_A_A, pi_A_B],  # Transitions from A
        [pi_B_A, pi_B_B]   # Transitions from B
    ])

    # Expected immediate reward vector R_pi[i] = E[r | s=i, pi]
    R_pi = np.array([
        pi_A_A * r_A_A + pi_A_B * r_A_B,  # Expected reward from A
        pi_B_A * r_B_A + pi_B_B * r_B_B   # Expected reward from B
    ])

    # Solve the system for V
    I = np.identity(2)
    A_matrix = I - gamma * P_pi
    V_pi = np.linalg.solve(A_matrix, R_pi)
    V_A, V_B = V_pi[0], V_pi[1]

    # 3. Calculate State-Action Values Q^pi
    # Q(s->s') = r(s->s') + gamma * V(s')
    Q_A_to_A = r_A_A + gamma * V_A
    Q_A_to_B = r_A_B + gamma * V_B

    # 4. Calculate Discounted State Occupancy d^pi
    # d^T = mu0^T * (I - gamma * P_pi)^-1
    A_inv = np.linalg.inv(A_matrix)
    d_pi = mu0 @ A_inv 
    d_A = d_pi[0]

    # 5. Compute the Gradients
    # dV/d(pi(a|s)) = d^pi(s) * Q^pi(s, a)
    grad_pi_A_A = d_A * Q_A_to_A
    grad_pi_A_B = d_A * Q_A_to_B

    # 6. Print the results using fractions
    f_d_A = Fraction(d_A).limit_denominator()
    f_Q_A_to_A = Fraction(Q_A_to_A).limit_denominator()
    f_Q_A_to_B = Fraction(Q_A_to_B).limit_denominator()
    f_grad_pi_A_A = Fraction(grad_pi_A_A).limit_denominator()
    f_grad_pi_A_B = Fraction(grad_pi_A_B).limit_denominator()

    print("The formula for the policy gradient is:")
    print("dV/d(pi(a|s)) = d^pi(s) * Q^pi(s, a)\n")
    
    print("Gradient for pi(A->A|A):")
    print(f"dV/d(pi(A->A|A)) = d^pi(A) * Q^pi(A, A->A)")
    print(f"= {f_d_A} * {f_Q_A_to_A}")
    print(f"= {f_grad_pi_A_A}\n")

    print("Gradient for pi(A->B|A):")
    print(f"dV/d(pi(A->B|A)) = d^pi(A) * Q^pi(A, A->B)")
    print(f"= {f_d_A} * {f_Q_A_to_B}")
    print(f"= {f_grad_pi_A_B}")
    
    # Store results for the final answer block
    global final_answer
    final_answer = f"The gradients are dV/d(pi(A->A|A)) = {f_grad_pi_A_A} and dV/d(pi(A->B|A)) = {f_grad_pi_A_B}"


solve_policy_gradient()
<<<The gradients are dV/d(pi(A->A|A)) = 162/169 and dV/d(pi(A->B|A)) = 54/169>>>