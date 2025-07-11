import numpy as np
from fractions import Fraction

def solve_policy_gradient():
    """
    Computes the policy gradients for the given MDP.
    """
    # Step 1: Define MDP and Policy Parameters
    print("Step 1: Define MDP and Policy Parameters")
    gamma = Fraction(1, 2)
    mu0 = np.array([Fraction(1, 4), Fraction(3, 4)], dtype=object)  # mu0(A), mu0(B)
    
    # Policy pi_1 probabilities
    pi_A_A = Fraction(1, 3)
    pi_A_B = Fraction(2, 3)
    pi_B_A = Fraction(1, 2)
    pi_B_B = Fraction(1, 2)
    
    # Reward values
    r_A_A = 1
    r_A_B = 0
    r_B_A = 1
    r_B_B = 0
    
    print(f"gamma = {gamma}")
    print(f"mu0 = [mu(A)={mu0[0]}, mu(B)={mu0[1]}]")
    print(f"pi(A->A|A) = {pi_A_A}, pi(A->B|A) = {pi_A_B}")
    print(f"pi(B->A|B) = {pi_B_A}, pi(B->B|B) = {pi_B_B}\n")

    # Step 2: Compute State-Value Function V^pi(s)
    print("Step 2: Compute State-Value Function V^pi(s)")
    
    # State transition matrix P_pi under the policy
    P_pi = np.array([
        [pi_A_A, pi_A_B],  # Transitions from state A
        [pi_B_A, pi_B_B]   # Transitions from state B
    ], dtype=object)

    # Expected immediate reward vector R_pi from each state
    R_pi = np.array([
        pi_A_A * r_A_A + pi_A_B * r_A_B,
        pi_B_A * r_B_A + pi_B_B * r_B_B
    ], dtype=object)

    # Solve the Bellman equation: V = (I - gamma * P_pi)^-1 * R_pi
    I = np.identity(2, dtype=object)
    mat_to_invert = I - gamma * P_pi
    
    # Manual inverse for 2x2 matrix of Fractions
    det = mat_to_invert[0, 0] * mat_to_invert[1, 1] - mat_to_invert[0, 1] * mat_to_invert[1, 0]
    inv_det = Fraction(1, 1) / det
    mat_inv = inv_det * np.array([
        [mat_to_invert[1, 1], -mat_to_invert[0, 1]],
        [-mat_to_invert[1, 0], mat_to_invert[0, 0]]
    ], dtype=object)

    V = mat_inv @ R_pi
    V_A, V_B = V[0], V[1]
    print(f"By solving the Bellman equations, we get:")
    print(f"V(A) = {V_A}")
    print(f"V(B) = {V_B}\n")

    # Step 3: Compute Action-Value Function Q^pi(s, a)
    print("Step 3: Compute Action-Value Function Q^pi(s, a)")
    Q_A_A = r_A_A + gamma * V_A
    Q_A_B = r_A_B + gamma * V_B
    print(f"Q(A, A->A) = r(A->A) + gamma * V(A) = {r_A_A} + {gamma} * {V_A} = {Q_A_A}")
    print(f"Q(A, A->B) = r(A->B) + gamma * V(B) = {r_A_B} + {gamma} * {V_B} = {Q_A_B}\n")

    # Step 4: Compute Discounted State Visitation Frequencies d^pi(s)
    print("Step 4: Compute Discounted State Visitation Frequencies d^pi(s)")
    d = mu0 @ mat_inv
    d_A, d_B = d[0], d[1]
    print("The frequencies d are calculated as d^T = mu0^T * (I - gamma * P_pi)^-1")
    print(f"d(A) = {d_A}")
    print(f"d(B) = {d_B}\n")

    # Step 5: Compute the Gradients
    print("Step 5: Compute the Gradients")

    # Gradient 1: dV / d(pi(A->A|A))
    # The derivative affects actions from state A. The formula simplifies to:
    # d(A) * [ (d/dp)pi(A->A|A)*Q(A,A->A) + (d/dp)pi(A->B|A)*Q(A,A->B) ]
    # Since pi(A->B|A) = 1 - pi(A->A|A), its derivative is -1.
    grad_A_A = d_A * (Q_A_A - Q_A_B)
    print(f"The first gradient is dV/d(pi(A->A|A)) = d(A) * [Q(A, A->A) - Q(A, A->B)]")
    print("Substituting the values:")
    print(f"= {d_A} * [{Q_A_A} - {Q_A_B}]")
    print(f"= {d_A} * [{Q_A_A - Q_A_B}]")
    print(f"= {grad_A_A}\n")
    
    # Gradient 2: dV / d(pi(A->B|A))
    # Similarly, the derivative of pi(A->A|A) w.r.t pi(A->B|A) is -1.
    grad_A_B = d_A * (Q_A_B - Q_A_A)
    print(f"The second gradient is dV/d(pi(A->B|A)) = d(A) * [Q(A, A->B) - Q(A, A->A)]")
    print("Substituting the values:")
    print(f"= {d_A} * [{Q_A_B} - {Q_A_A}]")
    print(f"= {d_A} * [{Q_A_B - Q_A_A}]")
    print(f"= {grad_A_B}\n")

    return grad_A_A, grad_A_B

if __name__ == '__main__':
    grad1, grad2 = solve_policy_gradient()
    # The final answer format
    print(f"<<<({grad1}, {grad2})>>>")
