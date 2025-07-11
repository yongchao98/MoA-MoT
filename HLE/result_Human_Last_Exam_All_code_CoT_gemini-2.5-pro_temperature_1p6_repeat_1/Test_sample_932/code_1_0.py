def analyze_collection_method():
    """
    Analyzes which insect tribes' immatures are unlikely to be collected by a beat-sheet method.

    The beat-sheet method is effective for insects on external foliage but not for
    internal borers, nest-dwellers, or those in soil/litter.
    """
    
    # Data: (index, tribe_name, immature_habitat, is_on_foliage)
    tribes_data = [
        (1, 'Apis', 'in a protected hive', False),
        (2, 'Melipotini', 'on foliage (caterpillars)', True),
        (3, 'Eupholini', 'inside wood/stems (borers)', False),
        (4, 'Acritini', 'in leaf litter/decaying matter', False),
        (5, 'Oxyptilini', 'on foliage (caterpillars)', True),
        (6, 'Dictyophorini', 'on stems/leaves (nymphs)', True),
        (7, 'Acanthocerini', 'in decaying wood/litter', False)
    ]

    unlikely_to_be_collected_indices = []

    print("Analyzing which tribes are unlikely to be collected via beat-sheet:")
    for index, name, habitat, is_on_foliage in tribes_data:
        if not is_on_foliage:
            unlikely_to_be_collected_indices.append(index)
            print(f"- Tribe {index} ({name}): Immatures are {habitat}. Unlikely to be collected.")
        else:
            print(f"- Tribe {index} ({name}): Immatures are {habitat}. Likely to be collected.")

    # Sort the indices in ascending order as requested
    unlikely_to_be_collected_indices.sort()

    # Format the final answer as a comma-separated string
    final_answer = ", ".join(map(str, unlikely_to_be_collected_indices))
    
    print("\nThe indices of the tribes whose immatures are unlikely to be collected are:")
    print(final_answer)

if __name__ == "__main__":
    analyze_collection_method()