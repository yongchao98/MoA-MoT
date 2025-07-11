import numpy as np

def print_state(name, state_vector):
    """Helper function to print a state vector nicely."""
    a = state_vector[0]
    b = state_vector[1]
    # Represent as a|0> + b|1>
    print(f"  {name} = ({a:.3f})|0> + ({b:.3f})|1>")

def print_equation(gate_matrix, initial_vector, result_vector, initial_name):
    """Helper function to print the matrix-vector multiplication."""
    print(f"Applying T gate to {initial_name}:")
    print("      [ {a:.2f}  {b:.2f} ]   [ {i0:.3f} ]   [ {r0:.3f} ]".format(
        a=gate_matrix[0, 0], b=gate_matrix[0, 1], i0=initial_vector[0], r0=result_vector[0]))
    print("      [ {c:.2f}  {d:.2f} ] @ [ {i1:.3f} ] = [ {r1:.3f} ]".format(
        c=gate_matrix[1, 0], d=gate_matrix[1, 1], i1=initial_vector[1], r1=result_vector[1]))
    print("-" * 30)

def main():
    """
    This script verifies that applying the T gate avoids a "death state"
    for all possible initial quantum lever states.
    """
    # Define constants
    SQRT2 = np.sqrt(2)

    # --- Define Basis States ---
    # Safe states in the computational basis
    s0 = np.array([1, 0], dtype=complex)
    s1 = np.array([0, 1], dtype=complex)

    # Other possible initial states
    s_minus = np.array([1/SQRT2, -1/SQRT2], dtype=complex)

    # Death states (and also possible initial states)
    s_i = np.array([1/SQRT2, 1j/SQRT2], dtype=complex)
    s_minus_i = np.array([1/SQRT2, -1j/SQRT2], dtype=complex)

    initial_states = {
        "|0>": s0,
        "|1>": s1,
        "|->": s_minus,
        "|i>": s_i,
        "|-i>": s_minus_i
    }

    # --- Define The Quantum Operation (T gate) ---
    # T gate = [[1, 0], [0, e^(i*pi/4)]]
    T_gate = np.array([[1, 0], [0, np.exp(1j * np.pi / 4)]], dtype=complex)

    print("Analyzing the T gate operation...\n")

    # --- Apply the gate and check the results ---
    for name, initial_state in initial_states.items():
        # Apply the gate: final_state = T * initial_state
        final_state = T_gate @ initial_state
        
        # Print the full operation
        print_equation(T_gate, initial_state, final_state, name)

        # Check if the final state is a death state.
        # This is true if the squared inner product |<death_state|final_state>|^2 is 1.
        prob_is_i = np.abs(np.vdot(s_i, final_state))**2
        prob_is_minus_i = np.abs(np.vdot(s_minus_i, final_state))**2

        print(f"Resulting state:")
        print_state(f"|ψ_final>", final_state)
        print(f"Probability of being |i>: {prob_is_i:.3f}")
        print(f"Probability of being |-i>: {prob_is_minus_i:.3f}")
        
        if np.isclose(prob_is_i, 1.0) or np.isclose(prob_is_minus_i, 1.0):
            print("Conclusion: DEATH STATE REACHED. This gate is not the solution.\n")
        else:
            print("Conclusion: SAFE. Not a death state.\n")

if __name__ == '__main__':
    main()