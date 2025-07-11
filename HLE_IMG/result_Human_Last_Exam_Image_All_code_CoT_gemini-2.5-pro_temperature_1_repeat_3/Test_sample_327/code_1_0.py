def solve_graphene_classification():
    """
    Analyzes and classifies three graphene nanoribbon band structures.
    """

    # The analysis for each plot is determined from the visual inspection of the band structures.
    # plot_data = [ {id, Edge, Width, Band} ]
    plots_data = [
        # Plot 1: Armchair, 8 dimer lines wide, Semiconducting
        {'id': 1, 'edge': 'A', 'width': 8, 'band': 1},
        # Plot 2: Armchair, 5 dimer lines wide, Metallic
        {'id': 2, 'edge': 'A', 'width': 5, 'band': 0},
        # Plot 3: Armchair, 6 dimer lines wide, Semiconducting
        {'id': 3, 'edge': 'A', 'width': 6, 'band': 1}
    ]

    classifications = []
    print("Analyzing each plot:")

    # Process each plot's data to form the classification string
    for plot in plots_data:
        # Create the classification string "Edge_Width_Band"
        classification = f"{plot['edge']}{plot['width']}{plot['band']}"
        classifications.append(classification)
        
        # Output the components for each part of the final string
        print(f"Plot {plot['id']}: Edge='{plot['edge']}', Width={plot['width']}, Band={plot['band']} => Classification='{classification}'")

    # Concatenate all individual classifications into the final result string
    final_string = "".join(classifications)

    print("\nFinal concatenated string:")
    print(final_string)

solve_graphene_classification()