import numpy as np

def run_simulation():
    """
    This script simulates the quantum trolley problem to find a safe operation.
    """
    # Define quantum states as complex vectors
    # |0> = [1, 0], |1> = [0, 1]
    psi_0 = np.array([1, 0], dtype=complex)
    psi_1 = np.array([0, 1], dtype=complex)
    # |-> = (|0> - |1>)/sqrt(2)
    psi_minus = (psi_0 - psi_1) / np.sqrt(2)
    # |i> = (|0> + i|1>)/sqrt(2) - Right track
    psi_i = (psi_0 + 1j * psi_1) / np.sqrt(2)
    # |-i> = (|0> - i|1>)/sqrt(2) - Left track
    psi_neg_i = (psi_0 - 1j * psi_1) / np.sqrt(2)

    # Dictionary of possible initial states with their names
    initial_states = {
        "|0>": psi_0,
        "|1>": psi_1,
        "|->": psi_minus,
        "|i>": psi_i,
        "|-i>": psi_neg_i,
    }

    # Define the T gate matrix: Rz(pi/4)
    # T = [[1, 0], [0, e^(i*pi/4)]]
    T_gate = np.array([[1, 0], [0, np.exp(1j * np.pi / 4)]], dtype=complex)
    
    print("Applying the T gate to all possible initial states:")
    print("-" * 50)

    all_safe = True
    # Loop through each possible initial state
    for name, psi_initial in initial_states.items():
        # Apply the T gate to the initial state
        psi_final = T_gate @ psi_initial

        # Calculate the probability of the final state being measured as |-i> (Left)
        # P(Left) = |<-i|psi_final>|^2
        prob_left = np.abs(np.vdot(psi_neg_i, psi_final))**2

        # Calculate the probability of the final state being measured as |i> (Right)
        # P(Right) = |<i|psi_final>|^2
        prob_right = np.abs(np.vdot(psi_i, psi_final))**2
        
        # Check if the outcome is deterministic (unsafe)
        # Using a tolerance for floating point comparisons
        if np.isclose(prob_left, 1.0) or np.isclose(prob_right, 1.0):
            all_safe = False

        print(f"Initial state: {name}")
        print(f"  P(Left) = |<-i|T|{name}>|^2 = {prob_left:.4f}")
        print(f"  P(Right) = |< i|T|{name}>|^2 = {prob_right:.4f}")
        print(f"  Outcome: {'UNSAFE' if not all_safe else 'Safe'}\n")

    print("-" * 50)
    if all_safe:
        print("Conclusion: The T gate is a safe operation for all possible initial states.")
    else:
        print("Conclusion: The T gate is not a safe operation.")

run_simulation()
<<<Y>>>