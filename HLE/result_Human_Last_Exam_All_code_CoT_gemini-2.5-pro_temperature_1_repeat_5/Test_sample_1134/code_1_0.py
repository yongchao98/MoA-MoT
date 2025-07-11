import numpy as np

def get_state_name(state_vector, state_name_map, tolerance=1e-5):
    """Finds the common name of a state vector if it exists in the map."""
    for name, vec in state_name_map.items():
        # Check for equality up to a global phase by calculating fidelity
        fidelity = np.abs(np.vdot(state_vector, vec))**2
        if np.isclose(fidelity, 1.0, atol=tolerance):
            return name
    # Fallback to numeric representation if no name is found
    return f"[{state_vector[0]:.2f}, {state_vector[1]:.2f}]"

def solve_trolley_problem():
    """
    Analyzes the quantum trolley problem to find the operation that avoids deaths.
    """
    # --- 1. Define States and Gates ---
    # Using the computational basis |0> = [1, 0], |1> = [0, 1]
    s_0 = np.array([1, 0], dtype=complex)
    s_1 = np.array([0, 1], dtype=complex)
    s_plus = (s_0 + s_1) / np.sqrt(2)
    s_minus = (s_0 - s_1) / np.sqrt(2)
    s_i = (s_0 + 1j * s_1) / np.sqrt(2)      # Right Track (Death State)
    s_neg_i = (s_0 - 1j * s_1) / np.sqrt(2)  # Left Track (Death State)

    # The chosen gate for the solution is Hadamard (H)
    H_gate = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)

    # --- 2. Set Up the Problem ---
    initial_states = {
        "|0⟩": s_0,
        "|1⟩": s_1,
        "|-⟩": s_minus,
        "|i⟩": s_i,
        "|-i⟩": s_neg_i,
    }

    death_states = {"|-i⟩": s_neg_i, "|i⟩": s_i}
    
    # A map to find pretty names for a given state vector
    all_named_states = {**initial_states, "|+⟩": s_plus}


    # --- 3. Run the Simulation ---
    print("Plan: Apply the Hadamard (H) gate to every possible initial state.")
    print("Goal: Ensure the final state is never |-i⟩ (Left Track) or |i⟩ (Right Track).\n")
    
    all_safe = True
    for name, initial_state in initial_states.items():
        # Apply the gate
        final_state = H_gate @ initial_state
        
        # Get symbolic names for printing the equation
        initial_name = get_state_name(initial_state, all_named_states)
        final_name = get_state_name(final_state, all_named_states)

        print(f"--- Scenario: Initial State is {initial_name} ---")
        # The problem asks to output the numbers in the final equation.
        # We show the equation symbolically for clarity.
        print(f"Equation: H {initial_name} = {final_name}")
        
        is_death_state = False
        for death_name, death_state_vec in death_states.items():
             # Check for physical equivalence (Fidelity == 1)
            fidelity = np.abs(np.vdot(final_state, death_state_vec))**2
            if np.isclose(fidelity, 1.0):
                print(f"RESULT: DANGER! Tram goes to track {death_name}.")
                is_death_state = True
                all_safe = False
                break
        
        if not is_death_state:
            print(f"RESULT: SAFE. The final state {final_name} is a superposition, not a deterministic track choice.")
        print("-" * 50)

    # --- 4. Conclusion ---
    print("\n--- FINAL CONCLUSION ---")
    if all_safe:
        print("The Hadamard (H) gate successfully avoids a deterministic track choice for all possible initial states.")
        print("This action places the tram in a superposition, ensuring no single track is chosen, thus saving everyone.")
        print("The correct choice is 'H'.")
    else:
        print("The chosen gate failed in at least one scenario.")

solve_trolley_problem()