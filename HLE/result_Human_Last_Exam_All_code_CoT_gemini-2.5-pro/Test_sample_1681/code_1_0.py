import numpy as np

def get_payoff_at_equilibrium():
    """
    Calculates the payoff for the Prisoner's Dilemma in a quantum setting.

    In the Eisert-Wilkens-Lewenstein (EWL) protocol for the Prisoner's Dilemma,
    with maximal entanglement (gamma = pi/2), a new Nash Equilibrium emerges
    where both players choose the 'Cooperate' (C) strategy.

    This function calculates the final payoff at this (C, C) equilibrium point.
    """

    # Payoff matrix values
    # R: Reward (C, C), S: Sucker (C, D), T: Temptation (D, C), P: Punishment (D, D)
    R = 5  # Payoff for Alice/Bob if (C, C)
    S = 0  # Payoff for Alice if (C, D)
    T = 7  # Payoff for Alice if (D, C)
    P = 1  # Payoff for Alice/Bob if (D, D)

    # Bob's payoffs for (C,D) and (D,C) are T and S respectively.
    # Alice's payoff is calculated as: P_A = R*p_cc + S*p_cd + T*p_dc + P*p_dd

    # In the (C, C) equilibrium, the quantum state collapses entirely to the
    # classical |CC> outcome. This is a known result for the EWL protocol
    # with maximal entanglement. We can demonstrate this by calculating the
    # probability amplitudes.

    # Player strategies are represented by unitary operators U(theta, phi).
    # Cooperate C is U(0, 0).
    # For Alice choosing C: theta_A = 0, phi_A = 0
    # For Bob choosing C: theta_B = 0, phi_B = 0
    theta_A, phi_A = 0, 0
    theta_B, phi_B = 0, 0
    
    c_A = np.cos(theta_A / 2)
    s_A = np.sin(theta_A / 2)
    c_B = np.cos(theta_B / 2)
    s_B = np.sin(theta_B / 2)

    # Amplitudes of the final state for maximal entanglement (gamma = pi/2)
    # c_ij corresponds to the amplitude of the state |ij> where i is Alice's
    # outcome and j is Bob's outcome (0 for C, 1 for D).
    c_00 = c_A * c_B * np.cos(phi_A + phi_B)
    c_01 = 1j * (c_B * s_A * np.cos(phi_B) - s_B * c_A * np.sin(phi_A))
    c_10 = 1j * (c_A * s_B * np.cos(phi_A) - s_A * c_B * np.sin(phi_B))
    c_11 = s_A * s_B + c_A * c_B * np.sin(phi_A + phi_B)

    # Probabilities are the squared magnitudes of the amplitudes
    p_cc = np.abs(c_00)**2
    p_cd = np.abs(c_01)**2
    p_dc = np.abs(c_10)**2
    p_dd = np.abs(c_11)**2
    
    # Calculate the payoff for Alice at the (C, C) equilibrium
    payoff_A = R * p_cc + S * p_cd + T * p_dc + P * p_dd
    
    # Calculate the payoff for Bob at the (C, C) equilibrium
    # Bob's Payoff: P_B = R*p_cc + T*p_cd + S*p_dc + P*p_dd
    payoff_B = R * p_cc + T * p_cd + S * p_dc + P * p_dd
    
    print("Quantum Prisoner's Dilemma Equilibrium:")
    print("-" * 40)
    print("The optimal strategy profile in the maximally entangled quantum game is (Cooperate, Cooperate).")
    print(f"Probabilities of outcomes at this equilibrium:")
    print(f"P(Cooperate, Cooperate) = {p_cc:.2f}")
    print(f"P(Cooperate, Defect)   = {p_cd:.2f}")
    print(f"P(Defect,   Cooperate) = {p_dc:.2f}")
    print(f"P(Defect,   Defect)   = {p_dd:.2f}")
    print("\nThis means both players are certain to cooperate, achieving the best mutual outcome.")
    
    print("\nThe payoff for each player is calculated as:")
    print(f"Payoff = ({R} * {p_cc:.2f}) + ({S} * {p_cd:.2f}) + ({T} * {p_dc:.2f}) + ({P} * {p_dd:.2f}) = {payoff_A:.2f}")
    print(f"The equilibrium point results in a payoff of ({int(payoff_A)}, {int(payoff_B)}) for the two players.")

get_payoff_at_equilibrium()