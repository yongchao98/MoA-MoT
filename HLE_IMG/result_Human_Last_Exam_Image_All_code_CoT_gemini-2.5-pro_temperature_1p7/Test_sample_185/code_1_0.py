import numpy as np
from fractions import Fraction

def solve_policy_gradient():
    """
    Computes the gradients of the value function with respect to policy parameters
    for the given MDP.
    """
    # 1. Define MDP parameters
    pi_AA = Fraction(1, 3)
    pi_AB = Fraction(2, 3)
    pi_BA = Fraction(1, 2)
    pi_BB = Fraction(1, 2)
    gamma = Fraction(1, 2)
    mu0_A = Fraction(1, 4)
    mu0_B = Fraction(3, 4)
    r_AA = 1
    r_AB = 0
    r_BA = 1
    r_BB = 0

    # 2. Set up and solve for V(A) and V(B)
    # Bellman equations in matrix form: (I - gamma * P_pi) * V = R_pi
    # (1 - gamma*pi_AA)V_A - gamma*pi_AB*V_B = pi_AA*r_AA + pi_AB*r_AB
    # -gamma*pi_BA*V_A + (1 - gamma*pi_BB)V_B = pi_BA*r_BA + pi_BB*r_BB

    A_matrix = np.array([
        [1 - gamma * pi_AA, -gamma * pi_AB],
        [-gamma * pi_BA, 1 - gamma * pi_BB]
    ])

    b_vector_V = np.array([
        pi_AA * r_AA + pi_AB * r_AB,
        pi_BA * r_BA + pi_BB * r_BB
    ])

    V = np.linalg.solve(A_matrix, b_vector_V)
    V_A, V_B = Fraction(V[0]).limit_denominator(), Fraction(V[1]).limit_denominator()

    # 3. Differentiate Bellman equations and solve for gradients wrt pi_AA
    # Let G_A = dV_A/dpi_AA and G_B = dV_B/dpi_AA
    # The constraint is pi_AB = 1 - pi_AA, so dpi_AB/dpi_AA = -1.
    # Differentiating the first equation wrt pi_AA:
    # (1-g*pi_AA)G_A - g*V_A - (g*(-1)*V_B + g*pi_AB*G_B) = r_AA
    # (1-g*pi_AA)G_A - g*pi_AB*G_B = r_AA + g*V_A - g*V_B
    #
    # Differentiating the second equation wrt pi_AA:
    # -g*pi_BA*G_A + (1-g*pi_BB)*G_B = 0
    
    b_vector_G = np.array([
        r_AA + gamma * V_A - gamma * V_B,
        0
    ])

    # The matrix of coefficients is the same A_matrix
    G = np.linalg.solve(A_matrix, b_vector_G)
    G_A, G_B = Fraction(G[0]).limit_denominator(), Fraction(G[1]).limit_denominator()

    # 4. Compute the final gradient for V(mu0) wrt pi_AA
    grad_V_pi_AA = mu0_A * G_A + mu0_B * G_B

    print("--- Gradient calculation for dV/d\u03C0(A\u2192A|A) ---")
    print(f"First, we compute the state values V(A) and V(B):")
    print(f"V(A) = {V_A}, V(B) = {V_B}")
    print("\nNext, we compute the gradients of the state values, dV(A)/d\u03C0(A\u2192A|A) and dV(B)/d\u03C0(A\u2192A|A):")
    print(f"dV(A)/d\u03C0(A\u2192A|A) = {G_A}")
    print(f"dV(B)/d\u03C0(A\u2192A|A) = {G_B}")
    print("\nFinally, we compute the gradient of the total value function:")
    print(f"dV(\u03BC\u2080)/d\u03C0(A\u2192A|A) = \u03BC\u2080(A) * dV(A)/d\u03C0(A\u2192A|A) + \u03BC\u2080(B) * dV(B)/d\u03C0(A\u2192A|A)")
    print(f"= {mu0_A} * {G_A} + {mu0_B} * {G_B}")
    print(f"= {mu0_A * G_A} + {mu0_B * G_B}")
    print(f"= {grad_V_pi_AA}")
    print("\n" + "="*30 + "\n")

    # 5. Compute the final gradient for V(mu0) wrt pi_AB
    # Since pi_AA + pi_AB = 1, d/dpi_AB = -d/dpi_AA
    grad_V_pi_AB = -grad_V_pi_AA
    
    print("--- Gradient calculation for dV/d\u03C0(A\u2192B|A) ---")
    print("Using the constraint \u03C0(A\u2192A|A) + \u03C0(A\u2192B|A) = 1, we apply the chain rule:")
    print("dV/d\u03C0(A\u2192B|A) = -dV/d\u03C0(A\u2192A|A)")
    print(f"= -({grad_V_pi_AA})")
    print(f"= {grad_V_pi_AB}")

    return grad_V_pi_AA, grad_V_pi_AB

if __name__ == '__main__':
    grad1, grad2 = solve_policy_gradient()
    # The final answer format is specified by the system.
    # The instruction requires <<<answer content>>> at the very end.
    # Let's prepare the string for it.
    final_answer = f"<<<{grad1}, {grad2}>>>"

solve_policy_gradient()