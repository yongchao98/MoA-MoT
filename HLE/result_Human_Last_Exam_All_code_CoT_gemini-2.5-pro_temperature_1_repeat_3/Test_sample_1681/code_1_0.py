import sympy as sp
from sympy.physics.quantum import TensorProduct

def solve_quantum_prisoners_dilemma():
    """
    Solves for the equilibrium payoff in the quantum Prisoner's Dilemma
    using the Eisert-Wilkens-Lewenstein (EWL) scheme.
    """
    # Define imaginary unit
    i = sp.I

    # 1. Define the payoff matrix from the problem description
    # (Cooperate, Cooperate), (Cooperate, Defect), (Defect, Cooperate), (Defect, Defect)
    # R is Reward, S is Sucker, T is Temptation, P is Punishment
    R, S, T, P = 5, 0, 7, 1

    # 2. Define the quantum mechanical operators for the EWL scheme
    
    # Basis vectors and Identity matrix
    C0 = sp.Matrix([1, 0])
    C1 = sp.Matrix([0, 1])
    I2 = sp.eye(2)
    
    # Pauli-X matrix
    sx = sp.Matrix([[0, 1], [1, 0]])

    # Entangling operator J and its conjugate transpose J_dag
    J = (1/sp.sqrt(2)) * (TensorProduct(I2, I2) + i * TensorProduct(sx, sx))
    J_dag = J.conjugate().transpose()

    # Initial state |00>
    psi_initial = TensorProduct(C0, C0)

    # Unitary operator for a player's strategy, as defined in the EWL paper
    def U_EWL(theta, phi):
        return sp.Matrix([
            [sp.exp(i*phi) * sp.cos(theta/2), sp.sin(theta/2)],
            [-sp.sin(theta/2), sp.exp(-i*phi) * sp.cos(theta/2)]
        ])

    # The 'miracle' quantum strategy Q is U(pi/2, pi/2)
    Q_op = U_EWL(sp.pi/2, sp.pi/2)
    
    # Assume both players adopt the quantum strategy Q
    U_A = Q_op
    U_B = Q_op

    # The combined strategy operator for both players
    U_full = TensorProduct(U_A, U_B)

    # 3. Calculate the final state of the system
    psi_final = J_dag * U_full * J * psi_initial

    # 4. Calculate the probabilities of the classical outcomes
    # The basis states correspond to (A's move, B's move)
    # |00> -> (Cooperate, Cooperate)
    # |01> -> (Cooperate, Defect)
    # |10> -> (Defect, Cooperate)
    # |11> -> (Defect, Defect)
    C00 = TensorProduct(C0, C0)
    C01 = TensorProduct(C0, C1)
    C10 = TensorProduct(C1, C1) # Typo in thought process, should be C1,C0
    C10 = TensorProduct(C1, C0)
    C11 = TensorProduct(C1, C1)
    
    p_CC = sp.Abs(psi_final.dot(C00.conjugate()))**2
    p_CD = sp.Abs(psi_final.dot(C01.conjugate()))**2
    p_DC = sp.Abs(psi_final.dot(C10.conjugate()))**2
    p_DD = sp.Abs(psi_final.dot(C11.conjugate()))**2
    
    # Simplify the probabilities to get clean numbers
    p_CC_s = sp.simplify(p_CC)
    p_CD_s = sp.simplify(p_CD)
    p_DC_s = sp.simplify(p_DC)
    p_DD_s = sp.simplify(p_DD)

    # 5. Calculate the expected payoff for Player A
    payoff_A = p_CC_s * R + p_CD_s * S + p_DC_s * T + p_DD_s * P
    
    # Print the equation for the final payoff
    print("The equilibrium is achieved when both players choose the quantum strategy Q.")
    print("This results in the following probabilities for the classical outcomes:")
    print(f"P(C,C) = {p_CC_s}")
    print(f"P(C,D) = {p_CD_s}")
    print(f"P(D,C) = {p_DC_s}")
    print(f"P(D,D) = {p_DD_s}\n")
    
    print("The payoff for each player is calculated as:")
    print(f"Payoff = P(C,C)*R + P(C,D)*S + P(D,C)*T + P(D,D)*P")
    print(f"Payoff = {p_CC_s}*{R} + {p_CD_s}*{S} + {p_DC_s}*{T} + {p_DD_s}*{P} = {sp.simplify(payoff_A)}")

solve_quantum_prisoners_dilemma()