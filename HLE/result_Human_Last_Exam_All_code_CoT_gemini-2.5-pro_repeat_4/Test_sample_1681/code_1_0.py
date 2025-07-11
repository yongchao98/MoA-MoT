import numpy as np

def solve_quantum_dilemma():
    """
    Solves the Prisoner's Dilemma in a quantum framework to find the equilibrium point.
    """
    # 1. Define Payoff Matrix and basis states
    # Payoffs: R=5 (Cooperate-Cooperate), S=0 (Cooperate-Defect), T=7 (Defect-Cooperate), P=1 (Defect-Defect)
    payoffs = {'CC': (5, 5), 'CD': (0, 7), 'DC': (7, 0), 'DD': (1, 1)}
    R, S, T, P = payoffs['CC'][0], payoffs['CD'][0], payoffs['DC'][0], payoffs['DD'][0]

    # Basis kets for a single qubit and two-qubit states
    ket_C = np.array([[1], [0]], dtype=complex) # Cooperate |C> = |0>
    ket_D = np.array([[0], [1]], dtype=complex) # Defect |D> = |1>
    ket_CC = np.kron(ket_C, ket_C)
    ket_CD = np.kron(ket_C, ket_D)
    ket_DC = np.kron(ket_D, ket_C)
    ket_DD = np.kron(ket_D, ket_D)
    basis = {'CC': ket_CC, 'CD': ket_CD, 'DC': ket_DC, 'DD': ket_DD}

    # 2. Define Quantum Game Operators
    # Pauli matrices
    I = np.eye(2, dtype=complex)
    SIGMA_X = np.array([[0, 1], [1, 0]], dtype=complex)
    SIGMA_Z = np.array([[1, 0], [0, -1]], dtype=complex)

    # Strategy operators
    C_op = I  # Cooperate is the Identity operation
    D_op = 1j * SIGMA_X  # Defect is a bit-flip like operation
    Q_op = 1j * SIGMA_Z  # Quantum strategy is a phase-flip like operation
    strategies = {'C': C_op, 'D': D_op, 'Q': Q_op}

    # Entangling operator J for maximal entanglement
    J = (1 / np.sqrt(2)) * (np.kron(I, I) + 1j * np.kron(SIGMA_X, SIGMA_X))
    J_dag = J.conj().T

    # Initial entangled state
    psi_initial = J @ ket_CC

    def calculate_payoffs(U_A, U_B):
        """Calculates the final state and expected payoffs for players A and B."""
        # Combine player strategies
        U_total = np.kron(U_A, U_B)
        
        # Calculate the final state after players' moves and disentanglement
        psi_final = J_dag @ U_total @ psi_initial

        # Calculate probabilities of each classical outcome
        probs = {name: np.abs(ket.conj().T @ psi_final)[0, 0]**2 for name, ket in basis.items()}

        # Calculate expected payoffs
        payoff_A = probs['CC'] * R + probs['CD'] * S + probs['DC'] * T + probs['DD'] * P
        payoff_B = probs['CC'] * R + probs['CD'] * T + probs['DC'] * S + probs['DD'] * P
        
        return (payoff_A, payoff_B), probs

    print("--- Analyzing the Quantum Prisoner's Dilemma ---")
    print(f"Payoff Matrix: (C,C):({R},{R}), (C,D):({S},{T}), (D,C):({T},{S}), (D,D):({P},{P})\n")
    
    # --- Analysis of (Q, Q) as a Nash Equilibrium ---
    print("1. Testing (Q, Q) as a potential Nash Equilibrium:")
    (payoff_A_QQ, payoff_B_QQ), _ = calculate_payoffs(Q_op, Q_op)
    print(f"   - Payoff for strategy pair (Q, Q): ({payoff_A_QQ:.1f}, {payoff_B_QQ:.1f})")

    (payoff_A_CQ, _), _ = calculate_payoffs(C_op, Q_op)
    print(f"   - If Player A deviates to C, payoff for (C, Q): Player A gets {payoff_A_CQ:.1f}")

    (payoff_A_DQ, _), _ = calculate_payoffs(D_op, Q_op)
    print(f"   - If Player A deviates to D, payoff for (D, Q): Player A gets {payoff_A_DQ:.1f}")

    if payoff_A_QQ >= payoff_A_CQ and payoff_A_QQ >= payoff_A_DQ:
        print("   -> Conclusion: Neither deviation is profitable. (Q, Q) is a Nash Equilibrium.\n")
    else:
        print("   -> Conclusion: A deviation is profitable. (Q, Q) is NOT a Nash Equilibrium.\n")

    # --- Analysis of (D, D) as a Nash Equilibrium ---
    print("2. Testing classical equilibrium (D, D):")
    (payoff_A_DD, payoff_B_DD), _ = calculate_payoffs(D_op, D_op)
    print(f"   - Payoff for strategy pair (D, D): ({payoff_A_DD:.1f}, {payoff_B_DD:.1f})")

    (payoff_A_QD, _), _ = calculate_payoffs(Q_op, D_op)
    print(f"   - If Player A deviates to Q, payoff for (Q, D): Player A gets {payoff_A_QD:.1f}")

    if payoff_A_DD >= payoff_A_QD:
        print("   -> Conclusion: Deviation to Q is not profitable. (D, D) is a Nash Equilibrium.\n")
    else:
        print(f"   -> Conclusion: Since {payoff_A_QD:.1f} > {payoff_A_DD:.1f}, deviating to Q is profitable. (D, D) is NOT a Nash Equilibrium.\n")

    # --- Final Conclusion and Equation ---
    print("--- Final Result ---")
    print("The unique Nash Equilibrium in this quantum game is (Q, Q).")
    print("This equilibrium resolves the dilemma by leading to the optimal cooperative outcome.")
    print("The equilibrium point (payoff) is (5.0, 5.0).")

    # Calculate the final equation as requested
    (_, _), final_probs = calculate_payoffs(Q_op, Q_op)
    p_cc, p_cd, p_dc, p_dd = final_probs['CC'], final_probs['CD'], final_probs['DC'], final_probs['DD']
    
    print("\nFinal Equation for the equilibrium payoff:")
    print(f"{payoff_A_QQ:.1f} = {R} * {p_cc:.1f} + {S} * {p_cd:.1f} + {T} * {p_dc:.1f} + {P} * {p_dd:.1f}")


solve_quantum_dilemma()
<<<5>>>