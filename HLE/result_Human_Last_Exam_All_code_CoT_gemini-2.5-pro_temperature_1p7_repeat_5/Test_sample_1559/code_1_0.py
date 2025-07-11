import math

def simulate_quantum_gates():
    """
    Simulates the passage of a classical bit through a sequence of quantum gates.
    """
    # The initial state is a classical bit 0.
    current_bit = 0
    num_cycles = 3

    print(f"Initial bit: {current_bit}\n")

    for i in range(1, num_cycles + 1):
        print(f"--- Cycle {i} ---")
        cycle_input_bit = current_bit

        # --- Gate A ---
        # (R1) Gate A puts its input into superposition of |0> and |1> states with equal probability.
        print(f"Input to Gate A: {cycle_input_bit}")
        # The actual state is now a superposition, which we will track conceptually.
        print("Output of Gate A: A quantum state in superposition.")

        # --- Gate B ---
        # (R2) Gate B performs a quantum measurement.
        # (R1) Gate A's output collapses to classical 1 if measured immediately afterwards.
        # The sequence is ABC, so B's measurement immediately follows A.
        print("Input to Gate B: The superposition state from Gate A.")
        bit_after_B = 1
        print(f"Output of Gate B: {bit_after_B} (due to decoherence from measurement as per rules R1 and R2).")

        # --- Gate C ---
        # (R3) Gate C applies a quantum-classical translation function.
        # Formula: |ψ⟩ → (|amplitude of |0⟩|² × 0 + |amplitude of |1⟩|² × 1)
        # We represent the classical input bit from Gate B as a quantum state to apply the formula.
        # Classical 0 is |0⟩ state -> amplitude of |0⟩ is 1.0, amplitude of |1⟩ is 0.0
        # Classical 1 is |1⟩ state -> amplitude of |0⟩ is 0.0, amplitude of |1⟩ is 1.0
        print(f"Input to Gate C: {bit_after_B}")
        
        amp_0 = 1.0 if bit_after_B == 0 else 0.0
        amp_1 = 1.0 if bit_after_B == 1 else 0.0

        # Apply the formula from R3
        bit_after_C = (amp_0**2 * 0) + (amp_1**2 * 1)

        # Output the calculation as requested
        print(f"Output of Gate C: The calculation is ({amp_0}**2 * 0) + ({amp_1}**2 * 1) = {int(bit_after_C)}")

        # Update the bit for the next cycle
        current_bit = int(bit_after_C)
        print("-" * 25)

    print(f"\nFinal classical output bit: {current_bit}")

simulate_quantum_gates()