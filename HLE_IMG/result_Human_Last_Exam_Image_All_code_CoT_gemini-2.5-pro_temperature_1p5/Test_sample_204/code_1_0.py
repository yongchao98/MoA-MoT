def sort_minerals_by_complexity():
    """
    Identifies minerals from images, classifies their silicate structure,
    and prints them in order of increasing structural complexity.
    """

    # Step 1: Define the specimens with their identified mineral type and silicate class.
    specimens = {
        'A': {'mineral': 'Pyroxene', 'structure': 'Inosilicate (Single Chain)'},
        'B': {'mineral': 'Mica/Chlorite', 'structure': 'Phyllosilicate (Sheet)'},
        'C': {'mineral': 'Olivine', 'structure': 'Nesosilicate (Island)'},
        'D': {'mineral': 'Quartz', 'structure': 'Tectosilicate (Framework)'}
    }

    # Step 2: Define the complexity hierarchy of silicate structures.
    # A lower number means a simpler structure.
    complexity_map = {
        'Nesosilicate (Island)': 1,
        'Inosilicate (Single Chain)': 2,
        'Phyllosilicate (Sheet)': 3,
        'Tectosilicate (Framework)': 4
    }

    print("Mineral Identification and Classification:")
    for key, data in specimens.items():
        print(f"Specimen {key} is {data['mineral']}, which is a {data['structure']}.")

    # Step 3: Sort the specimen keys based on the complexity of their silicate structure.
    sorted_keys = sorted(specimens.keys(), key=lambda k: complexity_map[specimens[k]['structure']])

    # Step 4: Print the final ordered list.
    print("\nOrder of Increasing Silicate Structure Complexity:")
    
    # Building the final equation string for display
    order_details = []
    for key in sorted_keys:
        mineral_info = f"{key} ({specimens[key]['structure']})"
        order_details.append(mineral_info)
        
    final_order_string = " < ".join(order_details)
    
    print(final_order_string)
    
    # Final answer as a simple comma-separated string for clarity.
    final_answer_simple = ", ".join(sorted_keys)
    print(f"\nFinal Answer: {final_answer_simple}")


sort_minerals_by_complexity()