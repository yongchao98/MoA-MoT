import numpy as np

def run_quantum_prisoners_dilemma():
    """
    Calculates and displays the equilibrium payoff for a quantum Prisoner's Dilemma.
    """
    # Payoff Matrix values
    # (C,C): (R,R) = (5,5)
    # (C,D): (S,T) = (0,7)
    # (D,C): (T,S) = (7,0)
    # (D,D): (P,P) = (1,1)
    R, S, T, P = 5, 0, 7, 1

    # Define basis states |C> = |0> and |D> = |1>
    C = np.array([1, 0])
    D = np.array([0, 1])

    # Basis states for the 2-qubit system
    CC = np.kron(C, C)
    CD = np.kron(C, D)
    DC = np.kron(D, C)
    DD = np.kron(D, D)

    # Define classical strategies as unitary operators
    U_C = np.identity(2)  # Cooperate = Identity
    U_D = np.array([[0, 1], [1, 0]])  # Defect = Pauli-X (sigma_x)

    # EWL framework: Create the entangling operator J for maximal entanglement (gamma = pi/2)
    I = np.identity(2)
    sigma_x = np.array([[0, 1], [1, 0]])
    J = (1 / np.sqrt(2)) * (np.kron(I, I) + 1j * np.kron(sigma_x, sigma_x))

    # Initial state is J|CC>
    psi_initial = J @ CC

    # A Nash Equilibrium for this payoff matrix is (Player 1: Defect, Player 2: Cooperate)
    # We create the corresponding operator U_D ⊗ U_C
    strategy_operator = np.kron(U_D, U_C)

    # Calculate the final state by applying the strategy operator
    psi_final = strategy_operator @ psi_initial

    # Calculate the probabilities of the four classical outcomes
    p_cc = np.abs(np.vdot(CC, psi_final))**2
    p_cd = np.abs(np.vdot(CD, psi_final))**2
    p_dc = np.abs(np.vdot(DC, psi_final))**2
    p_dd = np.abs(np.vdot(DD, psi_final))**2

    # Calculate the expected payoff for each player
    payoff_A = p_cc * R + p_cd * S + p_dc * T + p_dd * P
    payoff_B = p_cc * R + p_cd * T + p_dc * S + p_dd * P

    print("Chosen Equilibrium: (Player 1: Defect, Player 2: Cooperate)\n")
    print("Probabilities of classical outcomes:")
    print(f"P(C,C) = {p_cc:.2f}")
    print(f"P(C,D) = {p_cd:.2f}")
    print(f"P(D,C) = {p_dc:.2f}")
    print(f"P(D,D) = {p_dd:.2f}\n")

    print("Payoff Calculation:")
    print(f"Player 1 Payoff = P(C,C)*{R} + P(C,D)*{S} + P(D,C)*{T} + P(D,D)*{P}")
    print(f"                = {p_cc:.2f}*{R} + {p_cd:.2f}*{S} + {p_dc:.2f}*{T} + {p_dd:.2f}*{P} = {payoff_A:.2f}\n")

    print(f"Player 2 Payoff = P(C,C)*{R} + P(C,D)*{T} + P(D,C)*{S} + P(D,D)*{P}")
    print(f"                = {p_cc:.2f}*{R} + {p_cd:.2f}*{T} + {p_dc:.2f}*{S} + {p_dd:.2f}*{P} = {payoff_B:.2f}\n")
    
    print(f"The equilibrium point results in the payoff ({payoff_A:.2f}, {payoff_B:.2f}).")

run_quantum_prisoners_dilemma()
<<<3.5>>>