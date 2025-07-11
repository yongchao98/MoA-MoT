import numpy as np
from fractions import Fraction

def solve_mdp_gradient():
    """
    Solves for the policy gradients for the given MDP.
    The method involves setting up and solving two systems of linear equations:
    1. The Bellman equations to find the value function V.
    2. The differentiated Bellman equations to find the gradient of V.
    """
    # Step 1: Define MDP and Policy Parameters
    gamma = Fraction(1, 2)
    mu0 = {'A': Fraction(1, 4), 'B': Fraction(3, 4)}
    
    pi = {
        ('A', 'A'): Fraction(1, 3),
        ('A', 'B'): Fraction(2, 3),
        ('B', 'A'): Fraction(1, 2),
        ('B', 'B'): Fraction(1, 2)
    }
    
    r = {
        ('A', 'A'): 1,
        ('A', 'B'): 0,
        ('B', 'A'): 1,
        ('B', 'B'): 0
    }
    
    # Step 2: Calculate the state-value function V^pi(s)
    print("### Step 1: Calculate the state-value function V^pi ###")
    print("The Bellman equation V = R_pi + gamma * P_pi * V can be rewritten as (I - gamma * P_pi) * V = R_pi.")
    
    # Construct the policy transition matrix P_pi and expected reward vector R_pi
    P_pi = np.array([
        [pi[('A', 'A')], pi[('A', 'B')]],
        [pi[('B', 'A')], pi[('B', 'B')]]
    ], dtype=object)
    
    R_pi = np.array([
        pi[('A', 'A')] * r[('A', 'A')] + pi[('A', 'B')] * r[('A', 'B')],
        pi[('B', 'A')] * r[('B', 'A')] + pi[('B', 'B')] * r[('B', 'B')]
    ], dtype=object)

    # Solve the system (I - gamma * P_pi) * V = R_pi for V
    I = np.identity(2, dtype=object)
    M = I - gamma * P_pi
    
    # Use numpy's linalg.solve, which requires float type
    M_float = M.astype(float)
    R_pi_float = R_pi.astype(float)
    V_float = np.linalg.solve(M_float, R_pi_float)
    
    # Convert float results back to Fractions for precision
    V = np.array([Fraction(x).limit_denominator() for x in V_float])
    V_A, V_B = V[0], V[1]
    
    print(f"Solving for V gives: V(A) = {V_A}, V(B) = {V_B}\n")

    # Step 3: Calculate the derivative of V w.r.t p = pi(A->A|A)
    print("### Step 2: Calculate the gradient of the value function dV/dp ###")
    print("Using implicit differentiation, we solve (I - gamma * P_pi) * V' = R_pi' + gamma * P_pi' * V")
    
    # Derivative of R_pi wrt p
    R_pi_prime = np.array([Fraction(1), Fraction(0)], dtype=object)
    # Derivative of P_pi wrt p
    P_pi_prime = np.array([[Fraction(1), Fraction(-1)], [Fraction(0), Fraction(0)]], dtype=object)
    
    # Calculate the right-hand side vector C = R_pi' + gamma * P_pi' * V
    C = R_pi_prime + gamma * (P_pi_prime @ V)
    
    # Solve M * V_prime = C for V_prime
    C_float = C.astype(float)
    V_prime_float = np.linalg.solve(M_float, C_float)
    V_prime = np.array([Fraction(x).limit_denominator() for x in V_prime_float])
    V_A_prime, V_B_prime = V_prime[0], V_prime[1]

    print(f"Solving for V' = [dV(A)/dp, dV(B)/dp]^T gives: dV(A)/dp = {V_A_prime}, dV(B)/dp = {V_B_prime}\n")
    
    # Step 4: Compute the final gradient of V(mu0)
    print("### Step 3: Compute final gradients ###")
    
    # Gradient with respect to pi(A->A|A)
    grad_p_AA = mu0['A'] * V_A_prime + mu0['B'] * V_B_prime
    
    print("1. Gradient w.r.t. pi(A->A|A):")
    print(f"   d(V^pi(mu0))/d(pi(A->A|A)) = mu0(A) * dV(A)/dp + mu0(B) * dV(B)/dp")
    # Print the equation with numbers
    print(f"   = {mu0['A']} * {V_A_prime} + {mu0['B']} * {V_B_prime}")
    print(f"   = {mu0['A'] * V_A_prime} + {mu0['B'] * V_B_prime}")
    print(f"   = {grad_p_AA}\n")
    
    # Gradient with respect to pi(A->B|A)
    grad_p_AB = -grad_p_AA
    
    print("2. Gradient w.r.t. pi(A->B|A):")
    print(f"   Since pi(A->B|A) = 1 - pi(A->A|A), the gradient is the negative of the first.")
    print(f"   d(V^pi(mu0))/d(pi(A->B|A)) = {grad_p_AB}")
    
    final_answer = f"({grad_p_AA}, {grad_p_AB})"
    print(f"\n<<<{final_answer}>>>")
    
solve_mdp_gradient()