import sys
import io

# Redirect print to a string buffer to capture output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_coexpression_problem():
    """
    Analyzes plasmid combinations for co-expression suitability in E. coli.
    """
    # Define plasmid properties: origin of replication, resistance marker, and primary use.
    plasmids = {
        "pCDF-1b": {"origin": "CDF", "resistance": "spectinomycin", "type": "Expression"},
        "pET-28a(+)": {"origin": "ColE1", "resistance": "kanamycin", "type": "Expression"},
        "pGEX-T4-1": {"origin": "ColE1", "resistance": "ampicillin", "type": "Expression"},
        "pCDFDuet-1": {"origin": "CDF", "resistance": "spectinomycin", "type": "Duet Expression"},
        "pGEM-T": {"origin": "ColE1", "resistance": "ampicillin", "type": "Cloning (not for expression)"},
        "pET-15b": {"origin": "ColE1", "resistance": "ampicillin", "type": "Expression"},
        "pASK-IBA3": {"origin": "ColE1", "resistance": "chloramphenicol", "type": "Expression"},
        # Handle variants mentioned in incorrect options
        "pGEX-T4-1 (chloramphenicol)": {"origin": "ColE1", "resistance": "chloramphenicol", "type": "Expression"},
        "pET-28a(+) (ampicillin)": {"origin": "ColE1", "resistance": "ampicillin", "type": "Expression"},
    }

    # Define the answer choices
    choices = {
        "A": ["pCDF-1b", "pET-28a(+)"],
        "B": ["pET-28a(+)", "pGEX-T4-1"],
        "C": ["pCDFDuet-1", "pET-28a(+)"],
        "D": ["pET-28a(+)", "pGEX-T4-1 (chloramphenicol)"],
        "E": ["pGEM-T", "pCDF-1b"],
        "F": ["pET-15b", "pET-28a(+)"],
        "G": ["pASK-IBA3", "pET-28a(+) (ampicillin)"],
        "H": ["pCDFDuet-1"],
        "I": ["None"],
        "J": ["pGEX-T4-1", "pASK-IBA3"],
    }

    print("Analysis of Co-Expression Strategies:\n")
    print("The key principles for successful co-expression from two plasmids are:")
    print("1. Compatible Origins of Replication (e.g., ColE1 and CDF).")
    print("2. Different Antibiotic Resistance Markers.\n")

    valid_choices = []

    for choice, plasmid_names in choices.items():
        if choice == "I":
            continue

        print(f"--- Evaluating Option {choice}: {' & '.join(plasmid_names)} ---")
        
        if len(plasmid_names) == 1:
            p_name = plasmid_names[0]
            p_info = plasmids[p_name]
            if p_info["type"] == "Duet Expression":
                print(f"  Result: VALID. {p_name} is a single-plasmid Duet system designed for co-expression.")
                valid_choices.append(choice)
            else:
                print(f"  Result: INVALID. {p_name} alone is not a co-expression system.")
        else: # Two plasmids
            p1_name, p2_name = plasmid_names
            p1_info, p2_info = plasmids[p1_name], plasmids[p2_name]
            
            # Check 1: Origin Compatibility
            origins_compatible = p1_info["origin"] != p2_info["origin"]
            # Check 2: Resistance Marker Diversity
            markers_different = p1_info["resistance"] != p2_info["resistance"]
            # Check 3: Vector Functionality (no cloning vectors for expression)
            is_functional = "Cloning" not in p1_info["type"] and "Cloning" not in p2_info["type"]

            reasons = []
            if not origins_compatible:
                reasons.append(f"Incompatible Origins (both are {p1_info['origin']})")
            if not markers_different:
                reasons.append(f"Same Resistance Marker (both are {p1_info['resistance']})")
            if not is_functional:
                reasons.append(f"Improper Vector Type ({p1_name if 'Cloning' in p1_info['type'] else p2_name} is a cloning vector)")
            
            if origins_compatible and markers_different and is_functional:
                print(f"  Result: VALID.")
                print(f"  - Origins: {p1_info['origin']} and {p2_info['origin']} (Compatible)")
                print(f"  - Resistances: {p1_info['resistance']} and {p2_info['resistance']} (Different)")
                valid_choices.append(choice)
            else:
                print(f"  Result: INVALID. Reason(s): {', '.join(reasons)}.")
        print("")

    print("--- Conclusion ---")
    print(f"Valid options based on molecular biology principles are: {', '.join(valid_choices)}.")
    print("Let's determine the 'best' option:")
    print("- Option H (pCDFDuet-1) is an excellent single-plasmid solution.")
    print("- Option C (pCDFDuet-1 + pET-28a(+)) is a valid, but complex, system for co-expressing three proteins.")
    print("- Option A (pCDF-1b + pET-28a(+)) represents the most fundamental and widely used strategy for co-expression from two separate plasmids. It pairs a high-expression pET vector (ColE1 ori) with a compatible partner plasmid (CDF ori).")
    print("\nFor its robustness, common usage, and clear illustration of the plasmid compatibility principle, Option A is considered the best choice.")

# Execute the analysis and capture the output
solve_coexpression_problem()
# Get the captured output and print it to the actual console
output = captured_output.getvalue()
sys.stdout = old_stdout
print(output)
print("<<<A>>>")