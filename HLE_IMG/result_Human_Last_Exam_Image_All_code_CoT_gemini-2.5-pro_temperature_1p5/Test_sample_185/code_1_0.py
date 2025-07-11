import sympy

def solve_mdp_gradient():
    """
    Calculates the policy gradients for a given simple MDP.
    The code follows these steps:
    1.  Calculates the state-value functions V(s) by solving the Bellman equations.
    2.  Calculates the action-value functions Q(s, a).
    3.  Calculates the discounted state visitation frequencies d(s).
    4.  Computes the final gradients using the Policy Gradient Theorem.
    """
    # --- Define MDP parameters as SymPy Rational numbers for exact arithmetic ---
    gamma = sympy.Rational(1, 2)
    mu0_A = sympy.Rational(1, 4)
    mu0_B = sympy.Rational(3, 4)

    # Policy pi_1 parameters
    pi_A_A = sympy.Rational(1, 3)
    pi_A_B = sympy.Rational(2, 3)
    pi_B_A = sympy.Rational(1, 2)
    pi_B_B = sympy.Rational(1, 2)

    # Reward function r(s -> s')
    r_A_A = 1
    r_A_B = 0
    r_B_A = 1
    r_B_B = 0

    print("Step 1: Calculate the state value functions V(A) and V(B)")
    # V_A and V_B are the variables to solve for
    V_A, V_B = sympy.symbols('V_A V_B')

    # Bellman equations:
    # V(s) = sum_{a} pi(a|s) * [R(s,a) + gamma * sum_{s'} P(s'|s,a) * V(s')]
    eq1 = sympy.Eq(V_A, pi_A_A * (r_A_A + gamma * V_A) + pi_A_B * (r_A_B + gamma * V_B))
    eq2 = sympy.Eq(V_B, pi_B_A * (r_B_A + gamma * V_A) + pi_B_B * (r_B_B + gamma * V_B))
    
    solution_V = sympy.solve((eq1, eq2), (V_A, V_B))
    V_A_val = solution_V[V_A]
    V_B_val = solution_V[V_B]
    print(f"Solving the system of Bellman equations yields:")
    print(f"V(A) = {V_A_val}")
    print(f"V(B) = {V_B_val}\n")

    print("Step 2: Calculate the Q-values Q(s, a)")
    # Q(s, a) = R(s,a) + gamma * V(s')
    Q_A_A = r_A_A + gamma * V_A_val
    Q_A_B = r_A_B + gamma * V_B_val
    print(f"Q(A, A->A) = r(A->A) + gamma*V(A) = {r_A_A} + {gamma} * {V_A_val} = {Q_A_A}")
    print(f"Q(A, A->B) = r(A->B) + gamma*V(B) = {r_A_B} + {gamma} * {V_B_val} = {Q_A_B}\n")

    print("Step 3: Calculate the discounted state visitation frequencies d(s)")
    # Transition probabilities P(s'|s) under the policy
    P_A_A, P_B_A = pi_A_A, pi_A_B  # Transitions from A
    P_A_B, P_B_B = pi_B_A, pi_B_B  # Transitions from B
    
    # d_A and d_B are the variables to solve for
    d_A, d_B = sympy.symbols('d_A d_B')
    
    # Equations for state visitation frequencies: d(s) = mu(s) + gamma * sum_{s_prev} d(s_prev) * P(s|s_prev)
    eq_d1 = sympy.Eq(d_A, mu0_A + gamma * (d_A * P_A_A + d_B * P_A_B))
    eq_d2 = sympy.Eq(d_B, mu0_B + gamma * (d_A * P_B_A + d_B * P_B_B))
    
    solution_d = sympy.solve((eq_d1, eq_d2), (d_A, d_B))
    d_A_val = solution_d[d_A]
    d_B_val = solution_d[d_B]
    print(f"Solving the system of equations yields:")
    print(f"d(A) = {d_A_val}")
    print(f"d(B) = {d_B_val}\n")
    
    print("Step 4: Compute the policy gradients")

    # First gradient: derivative wrt pi(A->A|A)
    # The gradient of pi(A->A|A) is 1, and for pi(A->B|A) is -1. Others are 0.
    grad_1 = d_A_val * (1 * Q_A_A + (-1) * Q_A_B)
    
    print(f"--> Gradient dV/d(pi(A->A|A))")
    print(f"    = d(A) * [1 * Q(A, A->A) - 1 * Q(A, A->B)]")
    print(f"    = {d_A_val} * ({Q_A_A} - {Q_A_B})")
    print(f"    = {d_A_val} * ({Q_A_A - Q_A_B})")
    print(f"    = {grad_1}\n")
    
    # Second gradient: derivative wrt pi(A->B|A)
    # The gradient of pi(A->A|A) is -1, and for pi(A->B|A) is 1. Others are 0.
    grad_2 = d_A_val * ((-1) * Q_A_A + 1 * Q_A_B)

    print(f"--> Gradient dV/d(pi(A->B|A))")
    print(f"    = d(A) * [-1 * Q(A, A->A) + 1 * Q(A, A->B)]")
    print(f"    = {d_A_val} * ({Q_A_B} - {Q_A_A})")
    print(f"    = {d_A_val} * ({Q_A_B - Q_A_A})")
    print(f"    = {grad_2}\n")

    print("---------------- FINAL ANSWERS ----------------")
    print(f"The gradient dV/d(pi(A->A|A)) is: {grad_1}")
    print(f"The gradient dV/d(pi(A->B|A)) is: {grad_2}")

if __name__ == '__main__':
    solve_mdp_gradient()