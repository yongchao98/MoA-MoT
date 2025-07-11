import numpy as np

# Plan:
# 1. Define the payoff matrix values for the Prisoner's Dilemma.
# 2. Define the necessary quantum operators: Identity (I), Defect (X), entangling gate (J), and the quantum strategy (Q).
# 3. Set up the initial state |00>.
# 4. Simulate the game evolution for the strategy pair (Q, Q):
#    a. Create the entangled initial state: |psi_in> = J|00>
#    b. Apply players' strategies: |psi_f> = (Q ⊗ Q)|psi_in>
#    c. Disentangle the state for measurement: |psi_m> = J†|psi_f>
# 5. Calculate the probabilities of each classical outcome (CC, CD, DC, DD).
# 6. Calculate the final expected payoffs for both players and display the equation.

# 1. Payoff values
# (C,C) -> (R,R), (C,D) -> (S,T), (D,C) -> (T,S), (D,D) -> (P,P)
R = 5  # Reward
S = 0  # Sucker
T = 7  # Temptation
P = 1  # Punishment

print("Quantum Prisoner's Dilemma Equilibrium Calculation")
print("-" * 50)
print(f"Payoff Matrix Values: Reward(R)={R}, Sucker(S)={S}, Temptation(T)={T}, Punishment(P)={P}\n")

# 2. Quantum Operators
I = np.identity(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)  # Pauli-X, represents Defect (D)

# Entangling operator J for maximal entanglement (gamma = pi/2)
gamma = np.pi / 2
J = np.cos(gamma / 2) * np.kron(I, I) + 1j * np.sin(gamma / 2) * np.kron(X, X)

# The 'miracle' Quantum strategy Q
Q = np.array([[1j, 0], [0, -1j]], dtype=complex)

print("Players use the optimal quantum strategy, Q:\n", np.round(Q, 2))
print("\n--- Simulating the Game with Strategy (Q, Q) ---")

# 3. Initial State |00> = [1, 0, 0, 0]'
psi_0 = np.zeros(4, dtype=complex)
psi_0[0] = 1

# 4. Game Evolution
# a. Create entangled state
psi_in = J @ psi_0

# b. Apply players' strategies (Q ⊗ Q)
U_Alice = Q
U_Bob = Q
U_game = np.kron(U_Alice, U_Bob)
psi_f = U_game @ psi_in

# c. Disentangle the state
J_dagger = J.conj().T
psi_final = J_dagger @ psi_f

# 5. Calculate Probabilities
# Probabilities are the square of the absolute values of the amplitudes
probabilities = np.abs(psi_final)**2
P_CC = probabilities[0]  # Outcome |00> -> (Cooperate, Cooperate)
P_CD = probabilities[1]  # Outcome |01> -> (Cooperate, Defect)
P_DC = probabilities[2]  # Outcome |10> -> (Defect, Cooperate)
P_DD = probabilities[3]  # Outcome |11> -> (Defect, Defect)

print("\nFinal state after measurement corresponds to a probability distribution over classical outcomes:")
print(f"P(Cooperate, Cooperate) = {P_CC:.2f}")
print(f"P(Cooperate, Defect)   = {P_CD:.2f}")
print(f"P(Defect, Cooperate)   = {P_DC:.2f}")
print(f"P(Defect, Defect)     = {P_DD:.2f}")

# 6. Calculate Expected Payoffs
payoff_Alice = P_CC * R + P_CD * S + P_DC * T + P_DD * P
payoff_Bob = P_CC * R + P_CD * T + P_DC * S + P_DD * P

print("\n--- Calculating Final Payoffs ---")
print("The equilibrium point is the payoff vector (Payoff_Alice, Payoff_Bob).")

print("\nFinal Payoff Calculation for Alice:")
print(f"Payoff_Alice = P(CC)*R + P(CD)*S + P(DC)*T + P(DD)*P")
print(f"Payoff_Alice = ({P_CC:.2f} * {R}) + ({P_CD:.2f} * {S}) + ({P_DC:.2f} * {T}) + ({P_DD:.2f} * {P}) = {payoff_Alice:.2f}")

print("\nFinal Payoff Calculation for Bob:")
print(f"Payoff_Bob   = P(CC)*R + P(CD)*T + P(DC)*S + P(DD)*P")
print(f"Payoff_Bob   = ({P_CC:.2f} * {R}) + ({P_CD:.2f} * {T}) + ({P_DC:.2f} * {S}) + ({P_DD:.2f} * {P}) = {payoff_Bob:.2f}")

print("\nThus, the Nash Equilibrium in the quantum game results in a payoff of ({:.0f}, {:.0f}).".format(payoff_Alice, payoff_Bob))
print("This outcome is Pareto optimal and resolves the dilemma.")
print("<<<(5, 5)>>>")