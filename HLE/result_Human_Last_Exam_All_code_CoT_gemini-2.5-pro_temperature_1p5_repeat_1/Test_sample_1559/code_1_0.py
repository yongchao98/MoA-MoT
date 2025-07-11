import math

def simulate_quantum_system():
    """
    Solves the quantum-classical hybrid system problem by simulating the bit's state
    as it passes through the gate sequence ABCABCABC.
    """
    # The sequence of gates to be applied.
    sequence = "ABCABCABC"
    
    # Initial state is a classical bit 0.
    state = 0
    
    print(f"Starting simulation with initial classical bit: {state}")
    print(f"Applying gate sequence: {sequence}\n")
    
    # The sequence consists of 3 identical 'ABC' blocks.
    num_blocks = len(sequence) // 3
    
    for i in range(num_blocks):
        print(f"--- Processing block {i + 1} ('ABC') ---")
        print(f"Input to block: Classical Bit {state}")
        
        # --- Gate A ---
        # Rule R1: Gate A puts its input into superposition of |0> and |1> states
        # with equal probability, regardless of the input value.
        amp0 = 1 / math.sqrt(2)
        amp1 = 1 / math.sqrt(2)
        print(f"Step 1 (Gate A): State becomes a quantum superposition: {amp0:.4f}|0⟩ + {amp1:.4f}|1⟩.")

        # --- Gate B ---
        # Rule R2 defines Gate B as a measurement. Rule R1 specifies that when the
        # state from Gate A is measured, it collapses to a classical 1.
        state = 1
        print(f"Step 2 (Gate B): A measurement occurs, collapsing the state to Classical Bit {state}.")
        
        # --- Gate C ---
        # Rule R3 translates a state to a classical bit. If the input is already a
        # classical bit 'b', it is treated as the state |b⟩, and the formula
        # gives back 'b'.
        if state == 1:
            # For classical 1, the state is |1⟩, meaning 0|0⟩ + 1|1⟩
            c_amp0, c_amp1 = 0, 1
        else: # state == 0
            # For classical 0, the state is |0⟩, meaning 1|0⟩ + 0|1⟩
            c_amp0, c_amp1 = 1, 0

        output = (abs(c_amp0)**2 * 0) + (abs(c_amp1)**2 * 1)
        
        # This is the last step of the final block, so we show the detailed equation.
        if i == num_blocks - 1:
            print("Step 3 (Gate C): Applying the quantum-classical translation function.")
            print("   The final classical equation is: (|amplitude of |0⟩|² * 0) + (|amplitude of |1⟩|² * 1)")
            print(f"   Using amplitudes from state |{state}⟩ (amp|0⟩={c_amp0}, amp|1⟩={c_amp1}):")
            print(f"   Output = (|{c_amp0}|² * 0) + (|{c_amp1}|² * 1)")
            print(f"          = ({abs(c_amp0)**2} * 0) + ({abs(c_amp1)**2} * 1)")
            print(f"          = {abs(c_amp0)**2 * 0} + {abs(c_amp1)**2 * 1}")
            print(f"          = {int(output)}")
        else:
             print(f"Step 3 (Gate C): Processing classical bit {state}, output remains {int(output)}.")
        
        state = int(output)
        print("-" * 35)

    print("\nSimulation complete.")
    print(f"The final classical output bit is: {state}")

# Run the simulation
simulate_quantum_system()