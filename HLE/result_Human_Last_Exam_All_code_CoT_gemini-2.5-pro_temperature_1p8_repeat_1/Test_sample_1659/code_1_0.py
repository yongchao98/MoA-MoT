import cmath
import math
from collections import defaultdict

def apply_gate(state_vector, gate_fn):
    """Applies a given gate function to a state vector."""
    new_vector = defaultdict(complex)
    for basis_state, amplitude in state_vector.items():
        if amplitude == 0:
            continue
        # The gate function returns a list of (new_state, new_amplitude_factor)
        for new_state, amp_factor in gate_fn(basis_state):
            new_vector[new_state] += amplitude * amp_factor
    return new_vector

def h_gate(qubit_index):
    """Returns a function that applies the Hadamard gate to a specific qubit."""
    def gate_fn(basis_state):
        bit = basis_state[qubit_index]
        rest_of_state = list(basis_state)
        amp_factor = 1 / math.sqrt(2)
        
        if bit == '0':
            # |0> -> 1/sqrt(2) * (|0> + |1>)
            state0 = "".join(rest_of_state)
            rest_of_state[qubit_index] = '1'
            state1 = "".join(rest_of_state)
            return [(state0, amp_factor), (state1, amp_factor)]
        else: # bit == '1'
            # |1> -> 1/sqrt(2) * (|0> - |1>)
            state1 = "".join(rest_of_state)
            rest_of_state[qubit_index] = '0'
            state0 = "".join(rest_of_state)
            return [(state0, amp_factor), (state1, -amp_factor)]
    return gate_fn

def cnot_gate(control_index, target_index):
    """Returns a function that applies the CNOT gate."""
    def gate_fn(basis_state):
        if basis_state[control_index] == '1':
            new_state_list = list(basis_state)
            # Flip the target qubit
            new_state_list[target_index] = '1' if basis_state[target_index] == '0' else '0'
            return [("".join(new_state_list), 1)]
        else:
            return [(basis_state, 1)]
    return gate_fn

def ccnot_gate(control1_index, control2_index, target_index):
    """Returns a function that applies the Toffoli (CCNOT) gate."""
    def gate_fn(basis_state):
        if basis_state[control1_index] == '1' and basis_state[control2_index] == '1':
            new_state_list = list(basis_state)
            # Flip the target qubit
            new_state_list[target_index] = '1' if basis_state[target_index] == '0' else '0'
            return [("".join(new_state_list), 1)]
        else:
            return [(basis_state, 1)]
    return gate_fn

def main():
    """
    Simulates the quantum circuit and calculates the probability of measuring |100>.
    """
    # 1. Initial State: |psi_0> = |000>
    psi = {"000": complex(1.0, 0.0)}

    # 2. Apply a Hadamard gate to the first qubit (index 0)
    psi = apply_gate(psi, h_gate(0))
    
    # 3. Apply a CNOT gate (control=1st qubit, target=2nd qubit)
    psi = apply_gate(psi, cnot_gate(0, 1))

    # 4. Apply a Toffoli gate (controls=1st,2nd qubits, target=3rd qubit)
    psi = apply_gate(psi, ccnot_gate(0, 1, 2))

    # 5. Apply a second Hadamard gate to the first qubit
    psi = apply_gate(psi, h_gate(0))

    # Determine the probability of measuring |100>
    target_state = "100"
    amplitude = psi.get(target_state, 0)
    
    # Probability is the square of the magnitude of the amplitude
    probability = abs(amplitude)**2

    print(f"The final state |psi_4> is a superposition of states.")
    # The prompt asked for the numbers in the final equation.
    # The amplitude of |100> is part of this final state description.
    # The final probability calculation is the equation to be printed.
    print(f"The amplitude of the state |{target_state}> is {amplitude.real:.4f}")
    print(f"The probability of measuring |{target_state}> is |{amplitude.real:.4f}|^2 = {probability:.4f}")


if __name__ == "__main__":
    main()