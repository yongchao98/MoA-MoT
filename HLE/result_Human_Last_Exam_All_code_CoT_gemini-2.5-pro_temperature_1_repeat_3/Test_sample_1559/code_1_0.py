import math

def run_simulation():
    """
    Simulates the state of a bit passing through the ABCABCABC gate sequence.
    """
    # Initial state
    classical_bit = 0
    print(f"Starting simulation with initial classical bit: {classical_bit}\n")

    num_blocks = 3
    sequence = "ABC"

    for i in range(num_blocks):
        print(f"--- Processing Block {i + 1} ('{sequence}') ---")

        # --- Gate A ---
        # (R1) Gate A puts its input into superposition of |0> and |1> states with
        # equal probability but collapses to classical 1 if measured immediately afterwards.
        print(f"Step A: Input bit {classical_bit} enters Gate A.")
        # The specific quantum state is (1/sqrt(2))|0> + (1/sqrt(2))|1>
        print("         Gate A produces a superposition state: (1/√2)|0⟩ + (1/√2)|1⟩.")

        # --- Gate B ---
        # (R2) Gate B performs a quantum measurement that forces decoherence.
        # This measurement immediately follows Gate A, triggering the special condition in R1.
        print("Step B: The superposition state enters Gate B, which performs a measurement.")
        print("         Due to rule (R1), the measurement forces a deterministic collapse to 1.")
        classical_bit = 1
        print(f"         The bit is now classical: {classical_bit}.")

        # --- Gate C ---
        # (R3) Gate C applies a quantum-classical translation function.
        # We interpret the classical bit 1 as the quantum state |1>, which is 0|0> + 1|1>.
        print(f"Step C: Classical bit {classical_bit} enters Gate C.")
        print("         Rule (R3) is applied: output = (|amp of |0⟩|² × 0 + |amp of |1⟩|² × 1).")
        # For state |1>, the amplitude of |0> is 0 and the amplitude of |1> is 1.
        amp_0 = 0
        amp_1 = 1
        
        # Calculate the output as per the formula.
        output_value = (amp_0**2 * 0) + (amp_1**2 * 1)
        
        # Print the equation with the specific numbers used in the calculation.
        print(f"         Calculation: ({amp_0}**2 * 0) + ({amp_1}**2 * 1) = {int(output_value)}")
        
        classical_bit = int(output_value)
        print(f"         The output bit from Gate C is: {classical_bit}.")
        print(f"--- End of Block {i + 1}. Current bit is: {classical_bit} ---\n")

    print("="*30)
    print(f"Final classical output bit after the full sequence is: {classical_bit}")
    print("="*30)

if __name__ == "__main__":
    run_simulation()