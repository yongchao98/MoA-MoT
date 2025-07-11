import math

def calculate_quantum_prisoners_dilemma_equilibrium():
    """
    Calculates the equilibrium point for the Quantum Prisoner's Dilemma
    using the Eisert-Wilkens-Lewenstein (EWL) protocol.
    """
    # Payoff values from the problem description
    # (C,C) -> (R,R), (C,D) -> (S,T), (D,C) -> (T,S), (D,D) -> (P,P)
    R = 5  # Reward for mutual cooperation
    S = 0  # Sucker's payoff
    T = 7  # Temptation to defect
    P = 1  # Punishment for mutual defection

    # The quantum strategy Q is defined by the operator U(theta, phi)
    # The Nash Equilibrium in the quantum game is the strategy pair (Q, Q)
    # where Q = U(theta=pi/2, phi=pi/2).
    
    # Player A's strategy (Q)
    theta_A = math.pi / 2
    phi_A = math.pi / 2

    # Player B's strategy (Q)
    theta_B = math.pi / 2
    phi_B = math.pi / 2

    # Pre-calculate sine and cosine terms for the strategies
    # s_i = sin(theta_i / 2), c_i = cos(theta_i / 2)
    s_A = math.sin(theta_A / 2)
    c_A = math.cos(theta_A / 2)
    s_B = math.sin(theta_B / 2)
    c_B = math.cos(theta_B / 2)

    # Probabilities of the four classical outcomes (CC, CD, DC, DD)
    # These formulas are derived from the final state in the EWL protocol.
    # P_xy = |<xy|ψ_final>|^2
    
    # Probability of (Cooperate, Cooperate) outcome
    p_cc = (c_A * c_B * math.cos(phi_A + phi_B) - s_A * s_B)**2
    
    # Probability of (Cooperate, Defect) outcome for Player A
    p_cd = (c_A * s_B * math.cos(phi_A) - s_A * c_B * math.cos(phi_B))**2

    # Probability of (Defect, Cooperate) outcome for Player A
    p_dc = (c_A * s_B * math.sin(phi_A) - s_A * c_B * math.sin(phi_B))**2

    # Probability of (Defect, Defect) outcome
    p_dd = (s_A * s_B * math.cos(phi_A + phi_B) + c_A * c_B)**2
    
    # The sum of probabilities must be 1. Let's verify for the (Q,Q) case.
    # For (Q,Q): c_A=s_A=c_B=s_B=1/sqrt(2), phi_A=phi_B=pi/2
    # cos(pi)=-1, cos(pi/2)=0, sin(pi/2)=1
    # p_cc = (0.5*(-1) - 0.5)^2 = (-1)^2 = 1.0
    # p_cd = (0.5*0 - 0.5*0)^2 = 0.0
    # p_dc = (0.5*1 - 0.5*1)^2 = 0.0
    # p_dd = (0.5*(-1) + 0.5)^2 = 0.0
    # Sum = 1.0. The formulas are consistent for this case.

    # Calculate Player A's expected payoff at the (Q, Q) equilibrium
    payoff_A = R * p_cc + S * p_cd + T * p_dc + P * p_dd

    # Print the equation for the final payoff calculation
    print(f"{R} * {p_cc:.1f} + {S} * {p_cd:.1f} + {T} * {p_dc:.1f} + {P} * {p_dd:.1f} = {payoff_A:.1f}")

calculate_quantum_prisoners_dilemma_equilibrium()