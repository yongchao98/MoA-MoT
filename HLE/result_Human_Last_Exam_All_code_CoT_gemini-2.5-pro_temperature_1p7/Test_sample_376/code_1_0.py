def calculate_clonal_expansion():
    """
    Calculates the clonal expansion score based on CNVs, gene effects, and a repressor protein.
    """
    # Step 1: Define the primary tumor's genetic data
    cnvs = {
        'Chromosome 1': 3,  # Gain of 3
        'Chromosome 2': -2, # Loss of 2
        'Chromosome 3': 2   # Gain of 2
    }

    gene_data = {
        # Gene: [type, weight per copy change]
        'Oncogene A': ['oncogene', 0.5],
        'Tumor suppressor B': ['tsg', -0.7],
        'Oncogene C': ['oncogene', 0.4],
        'Tumor suppressor D': ['tsg', -0.6],
        'Oncogene E': ['oncogene', 0.3],
        'Tumor suppressor F': ['tsg', -0.5]
    }

    gene_locations = {
        'Chromosome 1': ['Oncogene A', 'Tumor suppressor D'],
        'Chromosome 2': ['Tumor suppressor B', 'Oncogene E'],
        'Chromosome 3': ['Oncogene C', 'Tumor suppressor F']
    }

    # Chromosomes where the tumor suppressor repressor is active
    repressor_active_on = {'Chromosome 1', 'Chromosome 2'}

    total_score = 0.0
    equation_components = []

    print("Calculating the clonal expansion score based on genetic events:\n")

    # Step 2 & 3: Iterate through chromosomes and apply rules
    for chromosome, copy_change in cnvs.items():
        for gene_name in gene_locations[chromosome]:
            gene_type, weight = gene_data[gene_name]

            contribution = 0.0

            # Oncogenes contribute on GAIN
            if gene_type == 'oncogene' and copy_change > 0:
                contribution = copy_change * weight
                print(f"Contribution from {gene_name} (gain): {copy_change} copies * {weight} = {contribution}")
            
            # Tumor suppressors (tsg) contribute on LOSS, unless repressed
            elif gene_type == 'tsg' and copy_change < 0:
                num_lost = abs(copy_change)
                if chromosome in repressor_active_on:
                    print(f"Contribution from {gene_name} (loss): Effect is NULLIFIED by repressor protein.")
                    contribution = 0.0
                else:
                    contribution = num_lost * weight
                    print(f"Contribution from {gene_name} (loss): {num_lost} copies * {weight} = {contribution}")
            
            if contribution != 0:
                total_score += contribution
                equation_components.append(str(contribution))

    # Step 4 & 5: Sum the scores and display the final result
    final_equation = " + ".join(equation_components)

    print("\n-----------------------------------------------------")
    print("Final Clonal Expansion Score Calculation:")
    print(f"Equation: {final_equation} = {total_score}")
    print("-----------------------------------------------------")


if __name__ == "__main__":
    calculate_clonal_expansion()
<<<2.3>>>