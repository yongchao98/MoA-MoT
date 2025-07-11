import numpy as np

# --- 1. Define Game Parameters ---
# Payoff matrix from the problem: [[(R,R), (S,T)], [(T,S), (P,P)]]
# R: Reward for mutual cooperation
# S: Sucker's payoff for cooperating while opponent defects
# T: Temptation to defect while opponent cooperates
# P: Punishment for mutual defection
R = 5
S = 0
T = 7
P = 1

# --- 2. Define Quantum States and Operators ---
# Computational basis vectors (kets)
ket_0 = np.array([[1], [0]], dtype=complex) # Represents Cooperate |C>
ket_1 = np.array([[0], [1]], dtype=complex) # Represents Defect |D>

# Basis states for the 2-qubit system using the tensor product (kron)
ket_CC = np.kron(ket_0, ket_0) # |CC> or |00>
ket_CD = np.kron(ket_0, ket_1) # |CD> or |01>
ket_DC = np.kron(ket_1, ket_0) # |DC> or |10>
ket_DD = np.kron(ket_1, ket_1) # |DD> or |11>

# Identity and Pauli matrices
I = np.eye(2, dtype=complex)
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

# Entangling operator J for creating a maximally entangled state
# J = exp(i * gamma/2 * sigma_x @ sigma_x) with gamma = pi/2
J = (1/np.sqrt(2)) * (np.kron(I, I) + 1j * np.kron(sigma_x, sigma_x))
# Disentangling operator is the Hermitian conjugate (dagger) of J
J_dag = J.conj().T

# Player strategies as unitary operators
C_op = I          # Cooperate strategy
D_op = 1j * sigma_y # Defect strategy (EWL protocol standard)
Q_op = 1j * sigma_z # Quantum "miracle" strategy

# --- 3. Set up and Run the Game for the (Q, Q) Equilibrium ---
# Initial state of the system is |CC>
psi_initial = ket_CC

# Create the maximally entangled state shared by the players
psi_entangled = J @ psi_initial

# Both players choose the quantum strategy Q
U_A = Q_op
U_B = Q_op

# The combined strategy operator is the tensor product of individual strategies
U_total = np.kron(U_A, U_B)

# Apply the strategies to the entangled state
psi_transformed = U_total @ psi_entangled

# Apply the disentangling operator to return to the computational basis before measurement
psi_final = J_dag @ psi_transformed

# --- 4. Calculate Probabilities of Outcomes ---
# The probability of a given outcome |xy> is the squared magnitude of the inner product: |<xy|psi_final>|^2
# Note: <xy| is the conjugate transpose of |xy>
P_CC = np.abs(ket_CC.conj().T @ psi_final)**2
P_CD = np.abs(ket_CD.conj().T @ psi_final)**2
P_DC = np.abs(ket_DC.conj().T @ psi_final)**2
P_DD = np.abs(ket_DD.conj().T @ psi_final)**2

# --- 5. Calculate and Print Final Payoffs ---
# Player A's expected payoff
payoff_A = P_CC * R + P_CD * S + P_DC * T + P_DD * P

# Player B's expected payoff (note swapped S and T)
payoff_B = P_CC * R + P_CD * T + P_DC * S + P_DD * P

print("Quantum Prisoner's Dilemma Equilibrium Point Calculation")
print("="*55)
print(f"Strategy Profile: (Player A: Quantum, Player B: Quantum)")
print(f"Resulting Probabilities:")
# .item() extracts the scalar value from the 1x1 numpy array
print(f"  P(Cooperate, Cooperate) = {P_CC.item():.2f}")
print(f"  P(Cooperate, Defect)    = {P_CD.item():.2f}")
print(f"  P(Defect,    Cooperate) = {P_DC.item():.2f}")
print(f"  P(Defect,    Defect)    = {P_DD.item():.2f}")
print("-"*55)
print("The equilibrium payoff for each player is calculated as:")
print("Payoff_A = P(CC)*R + P(CD)*S + P(DC)*T + P(DD)*P")
print("Payoff_B = P(CC)*R + P(CD)*T + P(DC)*S + P(DD)*P\n")

print("Player A Payoff Equation:")
print(f"Payoff = {P_CC.item():.2f} * {R} + {P_CD.item():.2f} * {S} + {P_DC.item():.2f} * {T} + {P_DD.item():.2f} * {P} = {payoff_A.item():.2f}")
print("\nPlayer B Payoff Equation:")
print(f"Payoff = {P_CC.item():.2f} * {R} + {P_CD.item():.2f} * {T} + {P_DC.item():.2f} * {S} + {P_DD.item():.2f} * {P} = {payoff_B.item():.2f}")
print(f"\nTherefore, the quantum equilibrium point is ({payoff_A.item():.0f}, {payoff_B.item():.0f}), which resolves the dilemma.")
