import numpy as np
from fractions import Fraction

def solve_mdp_gradients():
    """
    Solves for the policy gradients of the given MDP.
    """
    # Step 1: Define MDP parameters
    gamma = Fraction(1, 2)
    mu0 = np.array([Fraction(1, 4), Fraction(3, 4)]) # mu0(A), mu0(B)
    
    # Policy pi_1
    pi_AA = Fraction(1, 3) # pi(A->A|A)
    pi_AB = Fraction(2, 3) # pi(A->B|A)
    pi_BA = Fraction(1, 2) # pi(B->A|B)
    pi_BB = Fraction(1, 2) # pi(B->B|B)

    # Rewards r(s->s')
    r = {
        ('A', 'A'): Fraction(1), ('A', 'B'): Fraction(0),
        ('B', 'A'): Fraction(1), ('B', 'B'): Fraction(0)
    }

    # Step 2: Calculate V^pi(s)
    # The Bellman equation is V = R_pi + gamma * P_pi * V
    # which can be rewritten as (I - gamma * P_pi) * V = R_pi
    
    # Transition matrix P_pi
    P_pi = np.array([[pi_AA, pi_AB],
                     [pi_BA, pi_BB]])

    # Expected immediate reward vector R_pi
    R_pi_A = pi_AA * r[('A', 'A')] + pi_AB * r[('A', 'B')]
    R_pi_B = pi_BA * r[('B', 'A')] + pi_BB * r[('B', 'B')]
    R_pi = np.array([R_pi_A, R_pi_B])

    # Coefficient matrix for the linear system
    I = np.identity(2, dtype=object)
    A_mat = I - gamma * P_pi

    # Solve for V = [V(A), V(B)]
    V = np.linalg.solve(A_mat, R_pi)
    V_A, V_B = V[0], V[1]
    
    # --- Gradient w.r.t. pi(A->A|A) ---
    # Let theta = pi(A->A|A). Then pi(A->B|A) = 1 - theta.
    # Differentiating the Bellman equation (I - gamma*P)V = R w.r.t. theta:
    # (I - gamma*P)V' - gamma*(dP/dtheta)V = dR/dtheta
    # V' = (I - gamma*P)^-1 * (dR/dtheta + gamma*(dP/dtheta)V)

    # dP/dtheta
    dP_dtheta = np.array([[Fraction(1), Fraction(-1)],
                          [Fraction(0), Fraction(0)]])
    
    # dR/dtheta
    # R_A = theta*r(A->A) + (1-theta)*r(A->B) = theta*1 + (1-theta)*0 = theta
    # R_B is independent of theta.
    dR_dtheta = np.array([Fraction(1), Fraction(0)])

    # RHS of the V' equation
    rhs_V_prime = dR_dtheta + gamma * dP_dtheta @ V
    
    # Solve for V' = [dV(A)/dtheta, dV(B)/dtheta]
    V_prime = np.linalg.solve(A_mat, rhs_V_prime)
    dVAdpAA, dVBdpAA = V_prime[0], V_prime[1]
    
    # Final gradient calculation
    grad1 = mu0 @ V_prime

    # --- Print results for the first gradient ---
    print("--- Gradient calculation for ∂V(μ₀)/∂π(A→A|A) ---")
    print("\nStep 1: Calculate State Values V(A) and V(B)")
    print(f"Solving the Bellman equations, we get:")
    print(f"V(A) = {V_A.numerator}/{V_A.denominator}")
    print(f"V(B) = {V_B.numerator}/{V_B.denominator}")

    print("\nStep 2: Calculate Derivatives of State Values ∂V(s)/∂π(A→A|A)")
    print(f"Solving the differentiated Bellman equations, we get:")
    print(f"∂V(A)/∂π(A→A|A) = {dVAdpAA.numerator}/{dVAdpAA.denominator}")
    print(f"∂V(B)/∂π(A→A|A) = {dVBdpAA.numerator}/{dVBdpAA.denominator}")
    
    print("\nStep 3: Compute the final gradient")
    print("∂V(μ₀)/∂π(A→A|A) = μ₀(A) * ∂V(A)/∂π(A→A|A) + μ₀(B) * ∂V(B)/∂π(A→A|A)")
    
    term1_val = mu0[0] * dVAdpAA
    term2_val = mu0[1] * dVBdpAA
    print(f"= ({mu0[0].numerator}/{mu0[0].denominator}) * ({dVAdpAA.numerator}/{dVAdpAA.denominator}) + ({mu0[1].numerator}/{mu0[1].denominator}) * ({dVBdpAA.numerator}/{dVBdpAA.denominator})")
    print(f"= ({term1_val.numerator}/{term1_val.denominator}) + ({term2_val.numerator}/{term2_val.denominator})")
    print(f"= {grad1.numerator}/{grad1.denominator}")
    
    print("\n" + "="*50 + "\n")
    
    # --- Gradient w.r.t. pi(A->B|A) ---
    # Since pi(A->A|A) + pi(A->B|A) = 1, the derivative w.r.t. pi(A->B|A)
    # is the negative of the derivative w.r.t. pi(A->A|A).
    grad2 = -grad1
    dVAdpAB = -dVAdpAA
    dVBdpAB = -dVBdpAA

    # --- Print results for the second gradient ---
    print("--- Gradient calculation for ∂V(μ₀)/∂π(A→B|A) ---")
    print("\nSince π(A→A|A) + π(A→B|A) = 1, we have ∂/∂π(A→B|A) = -∂/∂π(A→A|A).")
    print(f"∂V(A)/∂π(A→B|A) = {dVAdpAB.numerator}/{dVAdpAB.denominator}")
    print(f"∂V(B)/∂π(A→B|A) = {dVBdpAB.numerator}/{dVBdpAB.denominator}")
    
    print("\n∂V(μ₀)/∂π(A→B|A) = μ₀(A) * ∂V(A)/∂π(A→B|A) + μ₀(B) * ∂V(B)/∂π(A→B|A)")
    
    term1_val_2 = mu0[0] * dVAdpAB
    term2_val_2 = mu0[1] * dVBdpAB
    print(f"= ({mu0[0].numerator}/{mu0[0].denominator}) * ({dVAdpAB.numerator}/{dVAdpAB.denominator}) + ({mu0[1].numerator}/{mu0[1].denominator}) * ({dVBdpAB.numerator}/{dVBdpAB.denominator})")
    print(f"= ({term1_val_2.numerator}/{term1_val_2.denominator}) + ({term2_val_2.numerator}/{term2_val_2.denominator})")
    print(f"= {grad2.numerator}/{grad2.denominator}")
    
    return grad1, grad2

if __name__ == '__main__':
    grad_AA, grad_AB = solve_mdp_gradients()
    # The final answer format
    print(f'<<<{grad_AA.numerator}/{grad_AA.denominator}, {grad_AB.numerator}/{grad_AB.denominator}>>>')
