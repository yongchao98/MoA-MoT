import numpy as np
from fractions import Fraction

def solve_policy_gradient():
    """
    Solves the policy gradient problem for the given MDP.
    """
    # --- 1. Define MDP Parameters ---
    gamma = 1/2
    mu0_A = 1/4
    mu0_B = 3/4
    pi_AA = 1/3  # pi(A->A | A)
    pi_AB = 2/3  # pi(A->B | A)
    pi_BA = 1/2  # pi(B->A | B)
    pi_BB = 1/2  # pi(B->B | B)
    r_AA = 1     # r(A->A)
    r_AB = 0     # r(A->B)
    r_BA = 1     # r(B->A)
    r_BB = 0     # r(B->B)

    # --- 2. Calculate State Values V^pi(s) ---
    # The Bellman equations form a system of linear equations:
    # V_A = pi_AA * (r_AA + gamma*V_A) + pi_AB * (r_AB + gamma*V_B)
    # => V_A = (1/3)*(1 + 1/2*V_A) + (2/3)*(0 + 1/2*V_B)
    # => (1 - gamma*pi_AA)*V_A - gamma*pi_AB*V_B = pi_AA*r_AA + pi_AB*r_AB
    # => (1 - 1/6)*V_A - (1/3)*V_B = 1/3 => (5/6)*V_A - (1/3)*V_B = 1/3 => 5*V_A - 2*V_B = 2
    #
    # V_B = pi_BA * (r_BA + gamma*V_A) + pi_BB * (r_BB + gamma*V_B)
    # => V_B = (1/2)*(1 + 1/2*V_A) + (1/2)*(0 + 1/2*V_B)
    # => -gamma*pi_BA*V_A + (1-gamma*pi_BB)*V_B = pi_BA*r_BA + pi_BB*r_BB
    # => (-1/4)*V_A + (1-1/4)*V_B = 1/2 => (-1/4)*V_A + (3/4)*V_B = 1/2 => -V_A + 3*V_B = 2
    
    coeffs_V = np.array([[5, -2], [-1, 3]])
    consts_V = np.array([2, 2])
    V_sol = np.linalg.solve(coeffs_V, consts_V)
    V_A, V_B = V_sol[0], V_sol[1]

    # --- 3. Calculate Action Values Q^pi(s, a) ---
    # Q(s, a') = r(s->a') + gamma * V(a')
    Q_AA = r_AA + gamma * V_A
    Q_AB = r_AB + gamma * V_B

    # --- 4. Calculate Discounted State Visitation Frequencies d^pi(s) ---
    # d(s) = mu0(s) + gamma * sum_s' d(s') * pi(s|s')
    # d_A = mu0_A + gamma * (d_A * pi_AA + d_B * pi_BA)
    # => (1 - gamma*pi_AA)*d_A - gamma*pi_BA*d_B = mu0_A
    # => (1 - 1/6)*d_A - (1/4)*d_B = 1/4 => (5/6)*d_A - (1/4)*d_B = 1/4 => 10*d_A - 3*d_B = 3
    #
    # d_B = mu0_B + gamma * (d_A * pi_AB + d_B * pi_BB)
    # => -gamma*pi_AB*d_A + (1-gamma*pi_BB)*d_B = mu0_B
    # => (-1/3)*d_A + (1-1/4)*d_B = 3/4 => (-1/3)*d_A + (3/4)*d_B = 3/4 => -4*d_A + 9*d_B = 9
    
    coeffs_d = np.array([[10, -3], [-4, 9]])
    consts_d = np.array([3, 9])
    d_sol = np.linalg.solve(coeffs_d, consts_d)
    d_A, d_B = d_sol[0], d_sol[1]

    # --- 5. Calculate Gradients ---
    # dV/dpi(a'|s) = d(s) * Q(s,a')
    grad_AA = d_A * Q_AA
    grad_AB = d_A * Q_AB

    # --- 6. Print Results ---
    print("This script computes the policy gradients using the Policy Gradient Theorem.")
    print("Formula: dV/dpi(a'|s) = d(s) * Q(s,a')\n")
    
    print("--- Intermediate Calculations ---")
    print(f"1. State Values: V(A) = {Fraction(V_A).limit_denominator()}, V(B) = {Fraction(V_B).limit_denominator()}")
    print(f"2. Action Values: Q(A,A) = {Fraction(Q_AA).limit_denominator()}, Q(A,B) = {Fraction(Q_AB).limit_denominator()}")
    print(f"3. State Visitation Frequency: d(A) = {Fraction(d_A).limit_denominator()}\n")

    print("--- Final Gradient Computations ---")
    # First Gradient
    print(f"Gradient 1: dV / d_pi(A->A|A) = d(A) * Q(A, A->A)")
    f_d_A = Fraction(d_A).limit_denominator()
    f_Q_AA = Fraction(Q_AA).limit_denominator()
    f_grad_AA = Fraction(grad_AA).limit_denominator()
    print(f"= {f_d_A} * {f_Q_AA} = {f_grad_AA}\n")
    
    # Second Gradient
    print(f"Gradient 2: dV / d_pi(A->B|A) = d(A) * Q(A, A->B)")
    f_Q_AB = Fraction(Q_AB).limit_denominator()
    f_grad_AB = Fraction(grad_AB).limit_denominator()
    print(f"= {f_d_A} * {f_Q_AB} = {f_grad_AB}")
    
    print("\n--- Final Answer ---")
    print(f"The gradient with respect to pi(A->A|A) is {f_grad_AA}")
    print(f"The gradient with respect to pi(A->B|A) is {f_grad_AB}")

solve_policy_gradient()