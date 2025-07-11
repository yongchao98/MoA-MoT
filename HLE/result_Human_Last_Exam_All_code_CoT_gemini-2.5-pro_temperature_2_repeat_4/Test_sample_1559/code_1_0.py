import math

def gate_A(classical_input):
    """
    (R1) Gate A puts its input into a superposition state |ψ⟩ with equal amplitudes
    for |0⟩ and |1⟩. The input value (0 or 1) doesn't alter the output state.
    """
    amplitude = 1 / math.sqrt(2)
    print("Gate A: The classical bit enters Gate A.")
    print(f"Gate A Output: A superposition state |ψ⟩ = {amplitude:.4f}|0⟩ + {amplitude:.4f}|1⟩.")
    return (amplitude, amplitude)

def gate_B(quantum_state_from_A):
    """
    (R2) Gate B performs a quantum measurement. According to (R1), a measurement
    immediately following Gate A forces a collapse to state |1⟩.
    """
    print("Gate B: A measurement is performed on the superposition state.")
    print("Gate B Output: Due to rule R1, the state deterministically collapses to |1⟩.")
    # The state |1⟩ is represented by amplitudes (0.0 for |0⟩, 1.0 for |1⟩).
    return (0.0, 1.0)

def gate_C(quantum_state_from_B, is_final_cycle):
    """
    (R3) Gate C applies a quantum-classical translation function to get a classical bit.
    Formula: (|amplitude of |0⟩|² × 0 + |amplitude of |1⟩|² × 1)
    """
    amp_0, amp_1 = quantum_state_from_B
    classical_output = (abs(amp_0)**2 * 0) + (abs(amp_1)**2 * 1)
    
    print("Gate C: The quantum state is converted back to a classical bit.")
    if is_final_cycle:
        # For the final step, output the equation with all its numbers.
        # We cast to int for cleaner presentation as the logic results in integers.
        print(f"Final Equation: (|{int(amp_0)}|² × 0 + |{int(amp_1)}|² × 1) = {int(classical_output)}")
    else:
        print(f"Calculation: ({amp_0:.1f}² × 0 + {amp_1:.1f}² × 1) = {int(classical_output)}")

    return int(classical_output)

# --- Simulation Start ---

# The initial classical bit is 0.
current_bit = 0
print(f"Initial classical bit: {current_bit}\n")

# The gate sequence is run three times.
num_cycles = 3
for i in range(num_cycles):
    print(f"--- Cycle {i + 1} of ABC ---")
    
    # Pass the current bit through the ABC sequence.
    state_after_A = gate_A(current_bit)
    state_after_B = gate_B(state_after_A)
    is_final = (i == num_cycles - 1)
    current_bit = gate_C(state_after_B, is_final_cycle=is_final)
    
    print(f"Output of cycle {i + 1}: Classical bit is now {current_bit}\n")

print(f"The final classical output bit after the full ABCABCABC sequence is: {current_bit}")