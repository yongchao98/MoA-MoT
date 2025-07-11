def analyze_antibody_response_mutants():
    """
    Analyzes mutant mouse groups to determine their expected impact on the
    production of high-affinity, somatically hypermutated antibodies.
    """
    groups = {
        'G1': {
            'gene': 'AID-(V18R)',
            'description': 'Activation-Induced Deaminase (AID) is the enzyme that directly mediates Somatic Hypermutation (SHM). A mutation in this essential enzyme will severely impair or abolish SHM and subsequent affinity maturation.',
            'affected': True
        },
        'G2': {
            'gene': 'CD40-KO',
            'description': 'The CD40 receptor on B cells is required to receive critical signals from T helper cells. Without this interaction, germinal centers do not form properly, preventing SHM and the generation of high-affinity antibodies.',
            'affected': True
        },
        'G3': {
            'gene': 'H2-IAd-(E137A/V142A)',
            'description': 'H2-IAd is the MHC Class II molecule that presents antigen peptides to CD4+ T helper cells. A mutation in the peptide-binding groove alters antigen presentation, thereby affecting the T cell help necessary for a robust antibody response.',
            'affected': True
        },
        'G4': {
            'gene': 'CD8-(V247D)',
            'description': 'CD8 is a co-receptor on cytotoxic T cells. The antibody response to a soluble protein like OVA depends on CD4+ T helper cells, not CD8+ T cells. Therefore, this mutation is not expected to significantly affect the outcome.',
            'affected': False
        },
        'G5': {
            'gene': 'H2-IAd-(T139A)',
            'description': 'Similar to G3, this mutation is in the MHC Class II molecule (H2-IAd). It will likely alter the presentation of OVA peptides to T helper cells, thus impacting the entire downstream process of affinity maturation.',
            'affected': True
        },
        'G6': {
            'gene': 'MyD88-KO',
            'description': 'MyD88 is an essential adaptor protein for the CpG adjuvant, which signals through TLR9. Without MyD88, the adjuvant effect is lost, leading to a significantly weaker immune response and lower antibody titers compared to wild-type.',
            'affected': True
        }
    }

    print("Step-by-step analysis of each mutant group:")
    print("-" * 75)

    affected_groups = []
    for group_id, data in groups.items():
        if data['affected']:
            affected_groups.append(group_id)
            status = "Affected"
            print(f"Group {group_id} ({data['gene']}): {status}")
            print(f"Reason: {data['description']}\n")

    print("-" * 75)
    print("Conclusion:")
    print("The mutant groups in which the titer of high-affinity OVA-specific antibodies")
    print("that have undergone SHM would be expected to be significantly different are:")

    # Printing each group ID as requested
    final_groups_str = ", ".join(affected_groups)
    print(final_groups_str)

# Execute the analysis
analyze_antibody_response_mutants()
<<<C>>>