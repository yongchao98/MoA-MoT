def analyze_collection_method():
    """
    Analyzes which insect tribes' immatures are unlikely to be collected
    using a beat-sheet method based on their life history.
    """
    # A dictionary storing the tribe name, the typical lifestyle of its immatures,
    # and a boolean indicating if they are likely to be collected by a beat-sheet.
    # True means likely, False means unlikely.
    tribe_data = {
        1: ("Apis", "Larvae live inside a protected hive/nest, fed by workers.", False),
        2: ("Melipotini", "Caterpillars (larvae) feed externally on plant leaves.", True),
        3: ("Eupholini", "Grubs (larvae) are typically internal borers in plant stems, roots, or wood.", False),
        4: ("Acritini", "Larvae are predators found in soil, leaf litter, or decaying organic matter.", False),
        5: ("Oxyptilini", "Caterpillars (larvae) feed externally on leaves, flowers, or buds.", True),
        6: ("Dictyophorini", "Nymphs are sap-feeders living on the exterior of plants.", True),
        7: ("Acanthocerini", "Larvae live in and feed on decaying wood or leaf litter.", False)
    }

    print("Analyzing the collection likelihood for each tribe's immatures:")
    print("-" * 80)

    unlikely_indices = []

    for index in sorted(tribe_data.keys()):
        tribe_name, lifestyle, is_collectible = tribe_data[index]
        if not is_collectible:
            unlikely_indices.append(index)
            status = "Unlikely"
        else:
            status = "Likely"
        
        print(f"Index {index} ({tribe_name}): {lifestyle} -> Collection Status: {status}")

    print("-" * 80)
    print("Based on the analysis, the immatures of the following tribes are unlikely to be collected:")
    
    # Format the final list of indices as a comma-separated string
    # The map(str, ...) converts each integer index to a string for joining.
    final_answer_string = ", ".join(map(str, unlikely_indices))
    
    print(final_answer_string)

# Run the analysis
analyze_collection_method()