import numpy as np

def solve_quantum_dilemma():
    """
    Calculates the equilibrium point for the Quantum Prisoner's Dilemma
    using the Eisert-Wilkens-Lewenstein (EWL) protocol with maximal entanglement.
    """
    # 1. Define the payoff values for Player 1 (Alice)
    # Payoffs: (CC, CD, DC, DD)
    payoff_cc = 5  # Reward for mutual cooperation
    payoff_cd = 0  # Sucker's payoff
    payoff_dc = 7  # Temptation to defect
    payoff_dd = 1  # Punishment for mutual defection

    print("Thinking Process:")
    print("1. The game is defined by the payoff matrix. Alice's payoffs are:")
    print(f"   (Cooperate, Cooperate): {payoff_cc}")
    print(f"   (Cooperate, Defect):   {payoff_cd}")
    print(f"   (Defect, Cooperate):   {payoff_dc}")
    print(f"   (Defect, Defect):     {payoff_dd}\n")

    # 2. Define the initial state with maximal entanglement (gamma = pi/2)
    # |psi_initial> = (1/sqrt(2)) * (|CC> + i|DD>)
    # In vector form |CC>=[1,0,0,0], |CD>=[0,1,0,0], |DC>=[0,0,1,0], |DD>=[0,0,0,1]
    psi_initial = (1 / np.sqrt(2)) * np.array([1, 0, 0, 1j], dtype=np.complex128)
    print("2. We set up a maximally entangled initial state for the players:")
    print(f"   |ψ_initial> = 1/√2 * (|CC> + i|DD>)\n")

    # 3. Define the optimal quantum strategy Q
    # Q = (i/sqrt(2)) * H, where H is the Hadamard matrix.
    H = np.array([[1, 1], [1, -1]], dtype=np.complex128)
    Q = (1j / np.sqrt(2)) * H
    print("3. Players can now choose a quantum strategy, Q:")
    print(f"   Q = {Q[0,0]:.2f} {Q[0,1]:.2f}")
    print(f"       {Q[1,0]:.2f} {Q[1,1]:.2f}\n")

    # 4. Calculate the final state when both players play Q
    # The combined operator is the tensor product U = Q ⊗ Q
    U = np.kron(Q, Q)
    # The final state is |ψ_final> = U * |ψ_initial>
    psi_final = U @ psi_initial
    print("4. Both players apply strategy Q. The final state |ψ_final> = (Q⊗Q)|ψ_initial> is calculated.\n")

    # 5. Calculate the probabilities of the four classical outcomes
    # Probability = |<outcome|ψ_final>|^2
    prob_cc = np.abs(psi_final[0])**2
    prob_cd = np.abs(psi_final[1])**2
    prob_dc = np.abs(psi_final[2])**2
    prob_dd = np.abs(psi_final[3])**2
    print("5. We find the probabilities of the classical outcomes from the final state:")
    print(f"   Prob(CC) = {prob_cc:.2f}")
    print(f"   Prob(CD) = {prob_cd:.2f}")
    print(f"   Prob(DC) = {prob_dc:.2f}")
    print(f"   Prob(DD) = {prob_dd:.2f}\n")

    # 6. Calculate the expected payoff for Alice at this (Q, Q) equilibrium
    expected_payoff = (payoff_cc * prob_cc +
                       payoff_cd * prob_cd +
                       payoff_dc * prob_dc +
                       payoff_dd * prob_dd)
    
    print("6. The equilibrium payoff is the expected value based on these probabilities.")
    print("\nFinal Calculation:")
    print("The payoff for a player at the quantum Nash Equilibrium (Q, Q) is calculated as follows:")
    
    # Final formatted output
    print(f"Payoff = P(CC)×Prob(CC) + P(CD)×Prob(CD) + P(DC)×Prob(DC) + P(DD)×Prob(DD)")
    print(f"Payoff = {payoff_cc} × {prob_cc:.2f} + {payoff_cd} × {prob_cd:.2f} + {payoff_dc} × {prob_dc:.2f} + {payoff_dd} × {prob_dd:.2f} = {expected_payoff:.2f}")

    return expected_payoff

# Run the calculation and print the final answer in the required format
result = solve_quantum_dilemma()
print(f"\n<<<Quantum Equilibrium Payoff: {result:.2f}>>>")
