import numpy as np

# Define the payoff matrix values
R = 5  # Reward for mutual cooperation
S = 0  # Sucker's payoff
T = 7  # Temptation to defect
P = 1  # Punishment for mutual defection

# --- 1. Define Quantum Operators ---
# 2x2 Identity matrix
I = np.identity(2, dtype=complex)
# Pauli matrices
SIGMA_X = np.array([[0, 1], [1, 0]], dtype=complex)
SIGMA_Z = np.array([[1, 0], [0, -1]], dtype=complex)

# Entangling operator J for maximal entanglement (gamma = pi/2)
J = (1 / np.sqrt(2)) * (I + 1j * np.kron(SIGMA_X, SIGMA_X))
# Adjoint (conjugate transpose) of J is the disentangling operator
J_dag = J.conj().T

# --- 2. Define Player Strategies ---
# In the quantum version, a new Nash Equilibrium strategy emerges.
# This quantum strategy, Q, resolves the dilemma.
Q = 1j * SIGMA_Z
# Both players will choose Q at equilibrium.
U_A = Q
U_B = Q
# The combined strategy operator is the tensor product of individual strategies.
U_combined = np.kron(U_A, U_B)

# --- 3. Simulate the Game ---
# Initial state |00> (Cooperate, Cooperate)
psi_0 = np.array([1, 0, 0, 0], dtype=complex)

# The game proceeds in three steps: entanglement, moves, disentanglement.
psi_final = J_dag @ U_combined @ J @ psi_0

# --- 4. Calculate Outcome Probabilities ---
# The probabilities are the squared magnitudes of the final state vector components.
# psi_final has components for |00>, |01>, |10>, |11>
prob_cc = np.abs(psi_final[0])**2
prob_cd = np.abs(psi_final[1])**2
prob_dc = np.abs(psi_final[2])**2
prob_dd = np.abs(psi_final[3])**2

# --- 5. Calculate Expected Payoffs ---
payoff_A = R * prob_cc + S * prob_cd + T * prob_dc + P * prob_dd
payoff_B = R * prob_cc + T * prob_cd + S * prob_dc + P * prob_dd

# --- 6. Print the Results ---
print("Quantum Prisoner's Dilemma Equilibrium Calculation\n")
print(f"Players' equilibrium strategy: Quantum strategy 'Q'")
print(f"Resulting outcome probabilities:")
print(f"  P(Cooperate, Cooperate) = {prob_cc:.1f}")
print(f"  P(Cooperate, Defect)   = {prob_cd:.1f}")
print(f"  P(Defect,   Cooperate) = {prob_dc:.1f}")
print(f"  P(Defect,   Defect)     = {prob_dd:.1f}\n")

print("Equilibrium Payoff Calculation (for each player):")
print(f"Payoff = R*P_CC + S*P_CD + T*P_DC + P*P_DD")
print(f"       = {R} * {prob_cc:.1f} + {S} * {prob_cd:.1f} + {T} * {prob_dc:.1f} + {P} * {prob_dd:.1f}")
print(f"       = {payoff_A:.1f}\n")
print(f"The equilibrium point is ({payoff_A:.1f}, {payoff_B:.1f}), which corresponds to the reward for mutual cooperation.")
