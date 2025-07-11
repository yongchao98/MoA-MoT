def find_unlikely_tribes():
    """
    Identifies insect tribes whose immatures are unlikely to be collected
    by a beat-sheet method based on their known habitats.
    """

    # Data: Tribe index, Tribe name, Immature habitat description.
    # Habitats are classified as either 'On Foliage' or 'Hidden/Internal'.
    tribes_data = {
        1: {'name': 'Apis', 'habitat': 'Hidden/Internal'},           # In hive cells
        2: {'name': 'Melipotini', 'habitat': 'On Foliage'},          # External caterpillar
        3: {'name': 'Eupholini', 'habitat': 'Hidden/Internal'},          # Internal borer/root feeder
        4: {'name': 'Acritini', 'habitat': 'Hidden/Internal'},           # In dung/carrion
        5: {'name': 'Oxyptilini', 'habitat': 'On Foliage'},          # External caterpillar
        6: {'name': 'Dictyophorini', 'habitat': 'On Foliage'},       # External nymph
        7: {'name': 'Acanthocerini', 'habitat': 'Hidden/Internal'}   # In soil/rotting wood
    }

    print("Analyzing which immatures are unlikely to be collected with a beat-sheet (which samples from foliage):")

    unlikely_indices = []
    for index, data in tribes_data.items():
        if data['habitat'] == 'Hidden/Internal':
            unlikely_indices.append(index)
            print(f"{index}) {data['name']}: Immatures are in a '{data['habitat']}' habitat, so they are UNLIKELY to be collected.")
        else:
            print(f"{index}) {data['name']}: Immatures are in an '{data['habitat']}' habitat, so they are LIKELY to be collected.")


    # Sort the indices for the final answer
    unlikely_indices.sort()

    # The final list is constructed from the individual numbers found
    print(f"\nThe indices of the unlikely tribes are {unlikely_indices[0]}, {unlikely_indices[1]}, {unlikely_indices[2]}, and {unlikely_indices[3]}.")
    
    # Format and print the final answer string as required.
    final_answer = ", ".join(map(str, unlikely_indices))
    print("\nFinal Answer:")
    print(final_answer)

find_unlikely_tribes()