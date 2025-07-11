import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a new string buffer
output_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = output_buffer


def solve_nmr_puzzle():
    """
    Analyzes 1H NMR data to identify the corresponding compound from a list of options.
    """
    # Step 1: Analyze the provided 1H NMR data.
    print("Step 1: Analyzing the 1H NMR Data")
    print("====================================")
    
    nmr_signals = {
        'aromatic': [(8.19, 1), (7.79, 1), (7.47, 1), (7.38, 1), (6.98, 1), (6.63, 1), (6.61, 1)],
        'aliphatic': [(4.19, 4), (3.63, 4), (3.21, 2), (2.83, 2), (1.98, 2)]
    }

    aromatic_protons = sum(integration for _, integration in nmr_signals['aromatic'])
    aliphatic_protons = sum(integration for _, integration in nmr_signals['aliphatic'])
    total_protons = aromatic_protons + aliphatic_protons
    
    print(f"The NMR spectrum shows {len(nmr_signals['aromatic'])} signals in the aromatic region.")
    print(f"Total integration of aromatic protons: {aromatic_protons}H")
    
    print(f"\nThe NMR spectrum shows {len(nmr_signals['aliphatic'])} signals in the aliphatic region.")
    print(f"Total integration of aliphatic protons: {aliphatic_protons}H")
    
    print("\nTotal observed protons in the spectrum:")
    # Print the equation as requested
    equation_parts = [str(integ) for group in nmr_signals.values() for _, integ in group]
    print(f"{' + '.join(equation_parts)} = {total_protons}H")
    print("\n")

    # Step 2: Analyze the candidate structures.
    print("Step 2: Analyzing the Chemical Structures")
    print("==========================================")
    
    # Analysis of Compound A
    print("--- Compound A ---")
    a_aromatic = 7  # 3 from quinoline + 4 from pyridine
    a_aliphatic = 14 # 6 from THQ aliphatic ring + 8 from piperazine
    a_nh = 1
    a_total = a_aromatic + a_aliphatic + a_nh
    print(f"Expected protons: {a_aromatic} (aromatic) + {a_aliphatic} (aliphatic) + {a_nh} (NH) = {a_total}H.")
    print("The single NH proton is acidic and often not observed in the spectrum due to exchange.")
    print(f"Expected observed protons (without NH): {a_total - a_nh}H.")
    print(f"Comparison: Expected {a_total - a_nh}H vs. Observed {total_protons}H -> MATCH.")
    print(f"Comparison (Aromatic): Expected {a_aromatic}H vs. Observed {aromatic_protons}H -> MATCH.")
    print(f"Comparison (Aliphatic): Expected {a_aliphatic}H vs. Observed {aliphatic_protons}H -> MATCH.")
    print("The number and integration of signals also match the structure of A.")
    
    # Analysis of other compounds
    print("\n--- Compounds B, D, E (Complexes) ---")
    print("These are [Zn(L)2]2+ complexes, containing two ligands. Their total proton count would be approximately double that of a single ligand (~44H), which is a significant mismatch with the observed 21H.")
    print("\n--- Compound C ---")
    c_aromatic = 8 # 3 from quinoline + 5 from phenyl
    c_aliphatic = 14
    c_total_no_nh = c_aromatic + c_aliphatic
    print(f"This has {c_aromatic} aromatic protons, which does not match the observed {aromatic_protons}H.")
    
    # Step 3: Conclusion.
    print("\nStep 3: Conclusion")
    print("====================")
    print("The NMR data is only consistent with the structure of Compound A.")
    print("The options provided are: A. B, B. D, C. E, D. C, E. A.")
    print("The correct choice is E, which corresponds to Compound A.")

solve_nmr_puzzle()

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = output_buffer.getvalue()

# Print the captured output to the actual console
print(output)
print("<<<E>>>")