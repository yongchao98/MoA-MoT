import numpy as np

def format_vector(v):
    """Formats a complex numpy vector for clean printing."""
    # Ensure it's a 1D array for consistent indexing
    v = np.ravel(v)
    parts = []
    for x in v:
        if np.isclose(x.imag, 0):
            parts.append(f"{x.real:6.3f}")
        else:
            parts.append(f"{x.real:6.3f}{x.imag:+.3f}j")
    return f"[{', '.join(parts)}]"

def run_simulation():
    """
    Simulates applying the Hadamard gate to the possible initial states
    of the quantum lever to determine if it's a safe operation.
    """
    # 1. Define the proposed quantum operation: the Hadamard gate (H)
    H_gate = 1/np.sqrt(2) * np.array([[1, 1], [1, -1]], dtype=complex)
    H_gate_str = "H"

    # 2. Define the possible initial states
    initial_states = {
        "|0>": np.array([[1], [0]], dtype=complex),
        "|1>": np.array([[0], [1]], dtype=complex),
        "|->": 1/np.sqrt(2) * np.array([[1], [-1]], dtype=complex),
        "|i>": 1/np.sqrt(2) * np.array([[1], [1j]], dtype=complex),
        "|-i>": 1/np.sqrt(2) * np.array([[1], [-1j]], dtype=complex)
    }

    # 3. Define the lethal output states
    lethal_i = 1/np.sqrt(2) * np.array([[1], [1j]], dtype=complex)
    lethal_neg_i = 1/np.sqrt(2) * np.array([[1], [-1j]], dtype=complex)

    print("Analyzing the effect of the Hadamard (H) gate on all possible initial states.")
    print("A 'Safe' result means the final state is not |-i> or |i>.")
    print("-" * 60)

    all_outcomes_are_safe = True
    for name, input_vec in initial_states.items():
        # Apply the gate to the state vector
        output_vec = H_gate @ input_vec

        # Check if the output is equivalent to a lethal state
        # We check if the inner product's magnitude squared is 1
        is_lethal_i = np.isclose(np.abs(np.vdot(lethal_i, output_vec))**2, 1)
        is_lethal_neg_i = np.isclose(np.abs(np.vdot(lethal_neg_i, output_vec))**2, 1)

        status = "Safe"
        if is_lethal_i or is_lethal_neg_i:
            status = "Lethal"
            all_outcomes_are_safe = False

        print(f"Operation: {H_gate_str} * {name}")
        print(f"Equation:")
        print(f"[[{H_gate[0,0].real:.3f}, {H_gate[0,1].real:.3f}]   {format_vector(input_vec)}   {format_vector(output_vec)}")
        print(f" [{H_gate[1,0].real:.3f}, {H_gate[1,1].real:.3f}]] *         =         ")
        print(f"Resulting state is {status}.\n")

    print("-" * 60)
    if all_outcomes_are_safe:
        print("Conclusion: The Hadamard (H) gate is the correct action. It guarantees safety in all scenarios.")
    else:
        print("Conclusion: The Hadamard (H) gate is not a safe choice.")

if __name__ == "__main__":
    run_simulation()