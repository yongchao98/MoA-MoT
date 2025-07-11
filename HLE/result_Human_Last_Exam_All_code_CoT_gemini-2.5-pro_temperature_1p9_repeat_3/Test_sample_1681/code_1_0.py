import numpy as np

def run_quantum_prisoners_dilemma():
    """
    Calculates and explains the equilibrium point for the quantum prisoner's dilemma.
    """
    # Payoff matrix from the problem description
    #               Player B: Cooperate (C)  Player B: Defect (D)
    # Player A: C | (5, 5)                  (0, 7)
    # Player A: D | (7, 0)                  (1, 1)
    payoff_matrix = {
        'CC': (5, 5), 'CD': (0, 7),
        'DC': (7, 0), 'DD': (1, 1)
    }

    # --- Step 1: Define Quantum Operators and States ---
    # Identity and Pauli matrices
    I = np.identity(2, dtype=complex)
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)

    # Player Strategies
    # Cooperate (C) is the identity operator
    U_C = I
    # Defect (D) is i * Pauli-X
    U_D = 1j * sigma_x

    # Initial state |00>
    psi_00 = np.array([1, 0, 0, 0], dtype=complex)

    # Entanglement operator J for maximal entanglement (gamma = pi/2)
    gamma = np.pi / 2
    J = np.cos(gamma / 2) * np.kron(I, I) + 1j * np.sin(gamma / 2) * np.kron(sigma_y, sigma_y)
    J_dag = J.conj().T

    # --- Step 2: Define a function to calculate payoffs ---
    def calculate_payoffs(U_A, U_B, name_A, name_B):
        # Create the combined strategy operator
        U_total = np.kron(U_A, U_B)
        
        # Calculate the final state vector
        psi_initial = J @ psi_00
        psi_mid = U_total @ psi_initial
        psi_final = J_dag @ psi_mid
        
        # Extract probabilities for |00>, |01>, |10>, |11>
        # Corresponds to (C,C), (C,D), (D,C), (D,D)
        probs = np.abs(psi_final)**2
        
        # Calculate expected payoffs
        payoff_A = (probs[0] * payoff_matrix['CC'][0] + 
                    probs[1] * payoff_matrix['CD'][0] +
                    probs[2] * payoff_matrix['DC'][0] + 
                    probs[3] * payoff_matrix['DD'][0])
        
        payoff_B = (probs[0] * payoff_matrix['CC'][1] + 
                    probs[1] * payoff_matrix['CD'][1] +
                    probs[2] * payoff_matrix['DC'][1] + 
                    probs[3] * payoff_matrix['DD'][1])

        print(f"Strategy: ({name_A}, {name_B})")
        print("Probabilities [P(CC), P(CD), P(DC), P(DD)]: "
              f"[{probs[0]:.2f}, {probs[1]:.2f}, {probs[2]:.2f}, {probs[3]:.2f}]")

        # Show the full calculation equation for Player A's payoff
        print(f"Player A Payoff Equation:\n"
              f"P_A = ({probs[0]:.2f} * {payoff_matrix['CC'][0]}) + "
              f"({probs[1]:.2f} * {payoff_matrix['CD'][0]}) + "
              f"({probs[2]:.2f} * {payoff_matrix['DC'][0]}) + "
              f"({probs[3]:.2f} * {payoff_matrix['DD'][0]}) = {payoff_A:.2f}")

        # Show the full calculation equation for Player B's payoff
        print(f"Player B Payoff Equation:\n"
              f"P_B = ({probs[0]:.2f} * {payoff_matrix['CC'][1]}) + "
              f"({probs[1]:.2f} * {payoff_matrix['CD'][1]}) + "
              f"({probs[2]:.2f} * {payoff_matrix['DC'][1]}) + "
              f"({probs[3]:.2f} * {payoff_matrix['DD'][1]}) = {payoff_B:.2f}")
        
        print(f"Resulting Payoffs (A, B): ({payoff_A:.2f}, {payoff_B:.2f})\n")
        
        return payoff_A, payoff_B

    # --- Step 3: Analyze the Equilibrium ---
    print("--- Finding the Quantum Nash Equilibrium ---\n")

    # Case 1: Both players Cooperate (C, C) - The potential equilibrium
    print("1. Calculate payoff for mutual cooperation (C, C):")
    eq_payoff_A, eq_payoff_B = calculate_payoffs(U_C, U_C, "Cooperate", "Cooperate")
    
    # Case 2: Player A deviates to Defect (D, C)
    print("2. Calculate payoff if Player A unilaterally defects (D, C):")
    dev_payoff_A, dev_payoff_B = calculate_payoffs(U_D, U_C, "Defect", "Cooperate")

    # --- Step 4: Conclusion ---
    print("--- Conclusion ---")
    print("In the quantum game, the payoff for mutual cooperation (C, C) is "
          f"({eq_payoff_A:.2f}, {eq_payoff_B:.2f}).")
    print("If Player A unilaterally deviates to Defect, their payoff becomes "
          f"{dev_payoff_A:.2f}.")
    if dev_payoff_A < eq_payoff_A:
        print(f"Since {dev_payoff_A:.2f} is less than {eq_payoff_A:.2f}, Player A has no "
              "incentive to defect. The same logic applies to Player B due to symmetry.")
        print("\nTherefore, (Cooperate, Cooperate) is a stable Nash Equilibrium, "
              "resolving the dilemma.")
        print("The equilibrium point is the payoff vector achieved at this equilibrium.")
    else:
        print("An error occurred, as deviation appears beneficial.")
    
    final_point = (round(eq_payoff_A), round(eq_payoff_B))
    print(f"\nFinal Equilibrium Point: {final_point}")


if __name__ == '__main__':
    run_quantum_prisoners_dilemma()
    # The direct answer format as requested. The equilibrium point is a vector in payoff space.
    print("<<<(5, 5)>>>")
