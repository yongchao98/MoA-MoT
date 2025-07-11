import numpy as np

def is_equivalent(state1, state2, tol=1e-9):
    """Checks if two states are the same up to a global phase."""
    # |<state2|state1>|^2 == 1 means they are equivalent
    return abs(np.vdot(state2, state1))**2 > 1 - tol

def run_trolley_simulation():
    """
    Simulates the quantum trolley problem to find a safe gate operation.
    """
    # Define basis states
    s0 = np.array([1, 0], dtype=complex)
    s1 = np.array([0, 1], dtype=complex)

    # Define the 5 possible initial states
    initial_states = {
        "|0>": s0,
        "|1>": s1,
        "|->": (s0 - s1) / np.sqrt(2),
        "|i>": (s0 + 1j * s1) / np.sqrt(2),
        "|-i>": (s0 - 1j * s1) / np.sqrt(2),
    }
    
    # Define the "kill" states
    kill_states = {
        "|i>": (s0 + 1j * s1) / np.sqrt(2),
        "|-i>": (s0 - 1j * s1) / np.sqrt(2),
    }

    # Define candidate quantum gates as matrices
    pi = np.pi
    gates = {
        "H": (1/np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex),
        "X": np.array([[0, 1], [1, 0]], dtype=complex),
        "S_dag": np.array([[1, 0], [0, -1j]], dtype=complex),
        "Z": np.array([[1, 0], [0, -1]], dtype=complex),
        "T": np.array([[1, 0], [0, np.exp(1j * pi / 4)]], dtype=complex),
        "T_dag": np.array([[1, 0], [0, np.exp(-1j * pi / 4)]], dtype=complex),
    }

    print("--- Starting Quantum Trolley Problem Simulation ---")
    
    solution_found = False
    for gate_name, gate_matrix in gates.items():
        print(f"\nTesting gate: {gate_name}")
        is_safe_gate = True
        
        for init_name, init_state in initial_states.items():
            final_state = gate_matrix @ init_state
            
            # Check if the final state is a kill state
            crashed = False
            for kill_name, kill_state in kill_states.items():
                if is_equivalent(final_state, kill_state):
                    print(f"  Applying {gate_name} to {init_name} results in {kill_name}. UNSAFE!")
                    
                    # Demonstrate the final state vector and kill state vector
                    # "output each number in the final equation!"
                    final_state_rounded = np.round(final_state, 3)
                    kill_state_rounded = np.round(kill_state, 3)
                    
                    # Normalize for comparison
                    phase = np.angle(final_state_rounded[np.nonzero(final_state_rounded)][0])
                    final_norm = np.round(final_state_rounded * np.exp(-1j*phase), 3)

                    phase_kill = np.angle(kill_state_rounded[np.nonzero(kill_state_rounded)][0])
                    kill_norm = np.round(kill_state_rounded * np.exp(-1j*phase_kill), 3)
                    
                    print(f"  Final state vector ~ {final_norm}")
                    print(f"  Kill state vector  ~ {kill_norm}")

                    is_safe_gate = False
                    crashed = True
                    break # No need to check other kill states
            if crashed:
                break # No need to check other initial states for this gate

        if is_safe_gate:
            print(f"  Gate {gate_name} is SAFE for all initial conditions.")
            print(f"SOLUTION: Applying the {gate_name} gate avoids all deaths.")
            solution_found = True
            
    if not solution_found:
        print("\nNo safe gate found in the tested list.")
        
    print("\n--- Simulation Complete ---")

# Run the simulation
run_trolley_simulation()