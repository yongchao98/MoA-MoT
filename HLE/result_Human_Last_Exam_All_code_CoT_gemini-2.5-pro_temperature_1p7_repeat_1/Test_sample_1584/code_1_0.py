import sys
# Redirect print to a string buffer to control the final output
from io import StringIO
old_stdout = sys.stdout
sys.stdout = mystdout = StringIO()


def analyze_coexpression_options():
    """
    Analyzes plasmid pairs for co-expression compatibility in E. coli.
    """
    # 1. Define a database of plasmid properties
    # Note: Some resistance markers in the options are intentionally incorrect
    # to test knowledge; we model those as separate entries.
    plasmid_db = {
        'pCDF-1b': {'origin': 'CloDF13', 'resistance': 'spectinomycin', 'notes': 'Standard vector.'},
        'pET-28a(+)': {'origin': 'ColE1', 'resistance': 'kanamycin', 'notes': 'High-level expression vector.'},
        'pGEX-T4-1': {'origin': 'ColE1', 'resistance': 'ampicillin', 'notes': 'Expression vector for GST-fusions.'},
        'pCDFDuet-1': {'origin': 'CloDF13', 'resistance': 'spectinomycin', 'notes': 'Specialized for co-expression (2 proteins).'},
        'pGEX-T4-1_alt': {'origin': 'ColE1', 'resistance': 'chloramphenicol'}, # As described in option D
        'pGEM-T': {'origin': 'ColE1', 'resistance': 'ampicillin', 'notes': 'Primarily a cloning vector, not for high-level expression.'},
        'pET-15b': {'origin': 'ColE1', 'resistance': 'ampicillin', 'notes': 'High-level expression vector.'},
        'pASK-IBA3': {'origin': 'ColE1', 'resistance': 'chloramphenicol', 'notes': 'Expression vector.'},
        'pET-28a(+)_alt': {'origin': 'ColE1', 'resistance': 'ampicillin'}, # As described in option G
    }

    # 2. Define the answer choices
    options = {
        'A': ['pCDF-1b', 'pET-28a(+)'],
        'B': ['pET-28a(+)', 'pGEX-T4-1'],
        'C': ['pCDFDuet-1', 'pET-28a(+)'],
        'D': ['pET-28a(+)', 'pGEX-T4-1_alt'],
        'E': ['pGEM-T', 'pCDF-1b'],
        'F': ['pET-15b', 'pET-28a(+)'],
        'G': ['pASK-IBA3', 'pET-28a(+)_alt'],
        'H': ['pCDFDuet-1'],
        'J': ['pGEX-T4-1', 'pASK-IBA3']
    }

    print("--- Analysis of Plasmid Co-expression Systems ---\n")
    print("Key requirements for a two-plasmid system:")
    print("1. Compatible Replication Origins (e.g., ColE1 and CloDF13)")
    print("2. Different Antibiotic Resistance Markers\n")
    
    viable_options = []

    for option, plasmids in options.items():
        print(f"--- Evaluating Option {option} ---")
        if len(plasmids) < 2:
            print("This option lists only one plasmid. A co-expression system with a separate protein of interest requires two.")
            print("Result: INCOMPLETE SYSTEM\n")
            continue

        p1_name, p2_name = plasmids
        p1 = plasmid_db[p1_name]
        p2 = plasmid_db[p2_name]

        print(f"Plasmid 1: {p1_name} (Origin: {p1['origin']}, Resistance: {p1['resistance']})")
        print(f"Plasmid 2: {p2_name} (Origin: {p2['origin']}, Resistance: {p2['resistance']})")

        # Check for compatibility
        origins_compatible = p1['origin'] != p2['origin']
        resistance_different = p1['resistance'] != p2['resistance']

        if not origins_compatible:
            print(f"Verdict: INCOMPATIBLE. Both plasmids have a {p1['origin']}-family origin.")
        elif not resistance_different:
            print(f"Verdict: INCOMPATIBLE. Both plasmids have the same resistance ({p1['resistance']}).")
        else:
            if 'cloning vector' in p1.get('notes', '') or 'cloning vector' in p2.get('notes', ''):
                print("Verdict: VIABLE, but suboptimal as it uses a cloning vector (e.g., pGEM-T) not designed for high expression.")
            else:
                print("Verdict: VIABLE. Origins are compatible and resistance markers are different.")
                viable_options.append(option)
        print("")

    # 3. Final recommendation
    print("--- Final Conclusion ---")
    if 'C' in viable_options and 'A' in viable_options:
        print("Both options A and C are viable systems based on compatibility rules.")
        print("However, Option C involves pCDFDuet-1, a vector specifically designed for co-expression. It contains two Multiple Cloning Sites (MCS), allowing expression of two different proteins (e.g., a chaperone complex like GroEL/GroES) from this single plasmid.")
        print("This specialized feature makes the system in Option C the most powerful and flexible choice.")
        print("\nTherefore, the best way to co-express a chaperone and protein of interest is represented by Option C.")
    elif 'C' in viable_options:
        print("Option C is the only fully viable and optimal system presented.")
    elif 'A' in viable_options:
        print("Option A is a viable system, but less flexible than one using a Duet vector.")
    else:
        print("None of the options describe a fully viable and optimal system for co-expression.")

# Run the analysis
analyze_coexpression_options()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(mystdout.getvalue())
print("<<<C>>>")