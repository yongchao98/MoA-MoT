import numpy as np

# Payoff matrix values
# (C,C) -> (R,R); (C,D) -> (S,T); (D,C) -> (T,S); (D,D) -> (P,P)
R = 5  # Reward for mutual cooperation
S = 0  # Sucker's payoff
T = 7  # Temptation to defect
P = 1  # Punishment for mutual defection

# Basis states |C> = |0> and |D_state> = |1>
C_state = np.array([1, 0])
D_state = np.array([0, 1])

# Two-qubit basis states
CC = np.kron(C_state, C_state)
CD = np.kron(C_state, D_state)
DC = np.kron(D_state, C_state)
DD = np.kron(D_state, D_state)

# Pauli Matrices
sigma_x = np.array([[0, 1], [1, 0]])
sigma_y = np.array([[0, -1j], [1j, 0]])
sigma_z = np.array([[1, 0], [0, -1]])
identity = np.identity(2)

# Entangling operator J for maximal entanglement (gamma = pi/2)
gamma = np.pi / 2
J = np.cos(gamma / 2) * np.identity(4) + 1j * np.sin(gamma / 2) * np.kron(sigma_x, sigma_x)
J_dagger = J.conj().T

# Player Strategies (Quantum Gates)
# C_op: Cooperate = Identity
C_op = identity
# D_op: Defect = U(pi, 0)
D_op = np.array([[0, 1], [-1, 0]])
# Q_op: Quantum = U(0, pi/2)
Q_op = np.array([[1j, 0], [0, -1j]])

# Initial entangled state
psi_in = J @ CC

def calculate_payoffs(U_A, U_B):
    """Calculates the final state and expected payoffs for players A and B."""
    # State after players apply their strategies
    psi_out = np.kron(U_A, U_B) @ psi_in
    
    # Final state after disentanglement
    psi_final = J_dagger @ psi_out
    
    # Probabilities of classical outcomes
    probs = {
        'CC': np.abs(CC.conj().T @ psi_final)**2,
        'CD': np.abs(CD.conj().T @ psi_final)**2,
        'DC': np.abs(DC.conj().T @ psi_final)**2,
        'DD': np.abs(DD.conj().T @ psi_final)**2
    }
    
    # Expected payoffs
    payoff_A = probs['CC'] * R + probs['CD'] * S + probs['DC'] * T + probs['DD'] * P
    payoff_B = probs['CC'] * R + probs['CD'] * T + probs['DC'] * S + probs['DD'] * P
    
    return payoff_A, payoff_B, probs

# --- Analysis ---
# 1. Payoff for (Q, Q) -> The ideal cooperative outcome
payoff_A_QQ, payoff_B_QQ, probs_QQ = calculate_payoffs(Q_op, Q_op)
print("Analyzing the Quantum Prisoner's Dilemma equilibrium:")
print("-" * 50)
print("Strategy (Alice, Bob) = (Quantum, Quantum)")
print(f"Probabilities: P(CC)={probs_QQ['CC']:.2f}, P(CD)={probs_QQ['CD']:.2f}, P(DC)={probs_QQ['DC']:.2f}, P(DD)={probs_QQ['DD']:.2f}")
print("Alice's Payoff Calculation:")
print(f"= {probs_QQ['CC']:.2f} * {R} (R) + {probs_QQ['CD']:.2f} * {S} (S) + {probs_QQ['DC']:.2f} * {T} (T) + {probs_QQ['DD']:.2f} * {P} (P) = {payoff_A_QQ:.2f}")
print("Bob's Payoff Calculation:")
print(f"= {probs_QQ['CC']:.2f} * {R} (R) + {probs_QQ['CD']:.2f} * {T} (T) + {probs_QQ['DC']:.2f} * {S} (S) + {probs_QQ['DD']:.2f} * {P} (P) = {payoff_B_QQ:.2f}")
print(f"Outcome Payoff: ({payoff_A_QQ:.2f}, {payoff_B_QQ:.2f})\n")

# 2. Check for unilateral deviation: Alice deviates from Q to D. Bob still plays Q.
payoff_A_DQ, payoff_B_DQ, probs_DQ = calculate_payoffs(D_op, Q_op)
print("Strategy (Alice, Bob) = (Defect, Quantum) -> Alice deviates")
print(f"Probabilities: P(CC)={probs_DQ['CC']:.2f}, P(CD)={probs_DQ['CD']:.2f}, P(DC)={probs_DQ['DC']:.2f}, P(DD)={probs_DQ['DD']:.2f}")
print("Alice's Payoff Calculation:")
print(f"= {probs_DQ['CC']:.2f} * {R} (R) + {probs_DQ['CD']:.2f} * {S} (S) + {probs_DQ['DC']:.2f} * {T} (T) + {probs_DQ['DD']:.2f} * {P} (P) = {payoff_A_DQ:.2f}")
print(f"Outcome Payoff: ({payoff_A_DQ:.2f}, {payoff_B_DQ:.2f})\n")

print("-" * 50)
print("Equilibrium Analysis:")
print(f"If both players choose Quantum (Q), their payoff is ({payoff_A_QQ:.2f}, {payoff_B_QQ:.2f}).")
print(f"However, if Alice unilaterally deviates to Defect (D), her payoff becomes {payoff_A_DQ:.2f}.")
print(f"Since {payoff_A_DQ:.2f} > {payoff_A_QQ:.2f}, Alice has an incentive to deviate.")
print("The same logic applies to Bob. Therefore, (Q,Q) is not a stable Nash Equilibrium.")
print("\nThe game between strategies Q and D becomes a new Prisoner's Dilemma.")
print("The only Nash Equilibrium is when both players choose to Defect (D,D).")
print("-" * 50)

# 3. Payoff for (D, D) -> The actual Nash Equilibrium of this quantum game
payoff_A_DD, payoff_B_DD, probs_DD = calculate_payoffs(D_op, D_op)
print("The Equilibrium Point: Strategy (Alice, Bob) = (Defect, Defect)")
print(f"Probabilities: P(CC)={probs_DD['CC']:.2f}, P(CD)={probs_DD['CD']:.2f}, P(DC)={probs_DD['DC']:.2f}, P(DD)={probs_DD['DD']:.2f}")
print("Equilibrium Payoff for Alice:")
print(f"= {probs_DD['CC']:.2f} * {R} + {probs_DD['CD']:.2f} * {S} + {probs_DD['DC']:.2f} * {T} + {probs_DD['DD']:.2f} * {P} = {payoff_A_DD:.2f}")
print("Equilibrium Payoff for Bob:")
print(f"= {probs_DD['CC']:.2f} * {R} + {probs_DD['CD']:.2f} * {T} + {probs_DD['DC']:.2f} * {S} + {probs_DD['DD']:.2f} * {P} = {payoff_B_DD:.2f}")
