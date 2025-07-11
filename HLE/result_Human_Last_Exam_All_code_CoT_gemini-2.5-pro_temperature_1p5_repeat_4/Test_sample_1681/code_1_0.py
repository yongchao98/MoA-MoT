import numpy as np

def print_equilibrium_calculation():
    """
    This function calculates and prints the payoff for the equilibrium point
    in the quantum prisoner's dilemma.
    """

    # Payoff matrix values: R(CC), S(CD), T(DC), P(DD)
    R, S, T, P = 5, 0, 7, 1
    payoffs_A = {'CC': R, 'CD': S, 'DC': T, 'DD': P}
    payoffs_B = {'CC': R, 'CD': T, 'DC': S, 'DD': P}

    # Define Pauli matrices and Identity
    ID = np.array([[1, 0], [0, 1]], dtype=complex)
    SIGMA_X = np.array([[0, 1], [1, 0]], dtype=complex)

    # EWL Protocol operators
    J = (1/np.sqrt(2)) * (np.kron(ID, ID) + 1j * np.kron(SIGMA_X, SIGMA_X))
    J_dagger = J.conj().T

    # Strategies
    C = ID  # Cooperate

    # The identified Nash Equilibrium in this quantum game is (Cooperate, Cooperate)
    U_A = C
    U_B = C

    # Initial state |00>
    psi_00 = np.array([1, 0, 0, 0], dtype=complex)

    # Calculate the final state after measurement
    # |psi_out> = J_dagger * (U_A kron U_B) * J * |00>
    psi_f = np.kron(U_A, U_B) @ J @ psi_00
    psi_out = J_dagger @ psi_f

    # Calculate probabilities of classical outcomes
    # P_xy = |<xy|psi_out>|^2
    prob_CC = np.abs(psi_out[0])**2
    prob_CD = np.abs(psi_out[1])**2
    prob_DC = np.abs(psi_out[2])**2
    prob_DD = np.abs(psi_out[3])**2

    # Calculate expected payoffs
    payoff_A = prob_CC * payoffs_A['CC'] + prob_CD * payoffs_A['CD'] + \
               prob_DC * payoffs_A['DC'] + prob_DD * payoffs_A['DD']
    payoff_B = prob_CC * payoffs_B['CC'] + prob_CD * payoffs_B['CD'] + \
               prob_DC * payoffs_B['DC'] + prob_DD * payoffs_B['DD']

    # Print the result as a detailed equation
    print("Equilibrium Point: (Cooperate, Cooperate)")
    print("-" * 40)
    print("Outcome Probabilities:")
    print(f"P(CC) = {prob_CC:.1f}, P(CD) = {prob_CD:.1f}, P(DC) = {prob_DC:.1f}, P(DD) = {prob_DD:.1f}")
    print("\nPayoff Calculation for Alice:")
    print(f"Payoff(A) = P(CC)*{payoffs_A['CC']} + P(CD)*{payoffs_A['CD']} + P(DC)*{payoffs_A['DC']} + P(DD)*{payoffs_A['DD']}")
    print(f"          = {prob_CC:.1f}*{payoffs_A['CC']} + {prob_CD:.1f}*{payoffs_A['CD']} + {prob_DC:.1f}*{payoffs_A['DC']} + {prob_DD:.1f}*{payoffs_A['DD']} = {payoff_A:.1f}")

    print("\nPayoff Calculation for Bob:")
    print(f"Payoff(B) = P(CC)*{payoffs_B['CC']} + P(CD)*{payoffs_B['CD']} + P(DC)*{payoffs_B['DC']} + P(DD)*{payoffs_B['DD']}")
    print(f"          = {prob_CC:.1f}*{payoffs_B['CC']} + {prob_CD:.1f}*{payoffs_B['CD']} + {prob_DC:.1f}*{payoffs_B['DC']} + {prob_DD:.1f}*{payoffs_B['DD']} = {payoff_B:.1f}")


print_equilibrium_calculation()