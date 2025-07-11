def analyze_antibody_response():
    """
    Analyzes the impact of genetic mutations on the production of high-affinity,
    somatically hypermutated (SHM) antibodies in mice.
    """
    # Dictionary mapping group ID to gene and biological rationale
    mouse_groups = {
        "G1": {
            "gene": "AID-(V18R)",
            "explanation": "AID is the enzyme essential for SHM. A mutation here prevents affinity maturation.",
            "is_affected": True
        },
        "G2": {
            "gene": "CD40-KO",
            "explanation": "CD40 is required for B cells to receive T-cell help, which is necessary for germinal center formation and SHM.",
            "is_affected": True
        },
        "G3": {
            "gene": "H2-IAd-(E137A/V142A)",
            "explanation": "H2-IAd is the MHC Class II molecule. This mutation impairs antigen presentation by B cells to T helper cells, disrupting the response.",
            "is_affected": True
        },
        "G4": {
            "gene": "CD8-(V247D)",
            "explanation": "CD8 is on cytotoxic T cells, not the T helper cells that assist B cells. This mutation does not directly affect the antibody response.",
            "is_affected": False
        },
        "G5": {
            "gene": "H2-IAd-(T139A)",
            "explanation": "Like G3, this MHC Class II mutation impairs the B cell's ability to get help from T cells.",
            "is_affected": True
        },
        "G6": {
            "gene": "MyD88-KO",
            "explanation": "MyD88 is required for the CpG adjuvant's effect via TLR9. Its absence leads to a significantly weaker immune response.",
            "is_affected": True
        }
    }

    print("Analysis of Mutant Groups for Altered Antibody Response:")
    print("-" * 60)

    affected_group_numbers = []

    for group_id, data in sorted(mouse_groups.items()):
        if data["is_affected"]:
            status = "AFFECTED"
            # Add the numeric part of the group ID to the list
            affected_group_numbers.append(group_id[1:])
        else:
            status = "UNAFFECTED"
        
        print(f"Group {group_id} ({data['gene']}): {status}")
        print(f"  Reason: {data['explanation']}\n")

    print("=" * 60)
    print("Conclusion: The groups with a significantly different antibody response are:")
    
    # Create the final output string, e.g., "G1, G2, G3, G5, G6"
    final_output = []
    for num in affected_group_numbers:
        final_output.append(f"G{num}")
    
    print(", ".join(final_output))


# Execute the analysis
analyze_antibody_response()