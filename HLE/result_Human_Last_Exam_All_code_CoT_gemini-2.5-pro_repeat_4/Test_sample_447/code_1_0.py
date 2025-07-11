import collections

def analyze_antibody_response_mutants():
    """
    Analyzes mutant mouse strains to predict their impact on the generation of
    high-affinity, somatically hypermutated (SHM) antibodies.
    """

    # Step 1: Define the function and essentiality of each gene/pathway for the
    # production of high-affinity, SHM-positive antibodies.
    # A value of 'True' for 'is_essential' means a KO/mutation is expected
    # to cause a significant difference compared to wild-type.
    protein_roles = {
        'AID': {
            'description': 'Activation-Induced Deaminase is the enzyme that directly mediates SHM. Its function is absolutely required.',
            'is_essential': True
        },
        'CD40': {
            'description': 'The CD40 receptor on B cells receives signals from T-helper cells. This interaction is critical for germinal center formation and affinity maturation.',
            'is_essential': True
        },
        'H2-IAd': {
            'description': 'This is an MHC Class II molecule that presents the OVA antigen to CD4+ T-helper cells. Proper presentation is required to initiate T-cell help.',
            'is_essential': True
        },
        'CD8': {
            'description': 'CD8 is a co-receptor on cytotoxic T cells, which are not directly involved in providing help to B cells for SHM in this context.',
            'is_essential': False
        },
        'MyD88': {
            'description': 'MyD88 is a key signaling adaptor for the CpG adjuvant (via TLR9). This adjuvant signal significantly boosts the B cell response and GC reaction.',
            'is_essential': True
        }
    }

    # Step 2: Define the mutant groups and the gene they affect.
    mutant_groups = collections.OrderedDict([
        ('G1', {'gene': 'AID', 'mutation': 'AID-(V18R)'}),
        ('G2', {'gene': 'CD40', 'mutation': 'CD40-KO'}),
        ('G3', {'gene': 'H2-IAd', 'mutation': 'H2-IAd-(E137A/V142A)'}),
        ('G4', {'gene': 'CD8', 'mutation': 'CD8-(V247D)'}),
        ('G5', {'gene': 'H2-IAd', 'mutation': 'H2-IAd-(T139A)'}),
        ('G6', {'gene': 'MyD88', 'mutation': 'MyD88-KO'})
    ])

    affected_groups = []

    print("Evaluating each mutant group:")
    print("-" * 60)

    # Step 3: Evaluate each group
    for group_id, info in mutant_groups.items():
        gene = info['gene']
        role = protein_roles[gene]
        
        print(f"Group {group_id} ({info['mutation']}):")
        print(f"  - Pathway Affected: {gene}")
        print(f"  - Rationale: {role['description']}")
        
        if role['is_essential']:
            print("  - Verdict: SIGNIFICANTLY DIFFERENT response expected.")
            affected_groups.append(group_id)
        else:
            print("  - Verdict: No significant difference expected.")
        print("-" * 60)

    # Step 4: Output the final results
    affected_groups.sort()
    affected_group_numbers = [int(g.replace('G', '')) for g in affected_groups]
    
    print("\n--- FINAL CONCLUSION ---")
    print("The groups expected to show a significantly different titer of high-affinity,")
    print("SHM-positive antibodies are those where an essential pathway is disrupted.")
    print(f"\nSelected Groups: {', '.join(affected_groups)}")

    print("\nThe final answer is composed of the numbers from the selected groups:")
    for num in affected_group_numbers:
        print(num)

# Execute the analysis
analyze_antibody_response_mutants()