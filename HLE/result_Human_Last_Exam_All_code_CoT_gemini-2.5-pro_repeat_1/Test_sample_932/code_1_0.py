def analyze_collection_method():
    """
    Analyzes which insect tribes' immatures are unlikely to be collected
    using a beat-sheet method based on their habitat.
    """

    # Data on the habitats of the immature stages for each tribe.
    # A beat-sheet primarily collects immatures living externally on foliage.
    tribes_data = {
        1: {'name': 'Apis', 'habitat': 'internal_hive', 'collectible': False},
        2: {'name': 'Melipotini', 'habitat': 'external_foliage', 'collectible': True},
        3: {'name': 'Eupholini', 'habitat': 'internal_borer', 'collectible': False},
        4: {'name': 'Acritini', 'habitat': 'detritus_dung', 'collectible': False},
        5: {'name': 'Oxyptilini', 'habitat': 'external_foliage', 'collectible': True},
        6: {'name': 'Dictyophorini', 'habitat': 'external_foliage', 'collectible': True},
        7: {'name': 'Acanthocerini', 'habitat': 'soil_wood', 'collectible': False}
    }

    print("The beat-sheet method collects insects living on the exterior of plants.")
    print("We must find the tribes whose immatures are not found on foliage.")
    print("-" * 30)

    unlikely_indices = []
    for index, data in tribes_data.items():
        if not data['collectible']:
            unlikely_indices.append(index)
            print(f"{index}) {data['name']}: Immatures are in habitat '{data['habitat']}'. Unlikely to be collected.")
        else:
            print(f"{index}) {data['name']}: Immatures are in habitat '{data['habitat']}'. Likely to be collected.")

    # The final answer requires the indices sorted in ascending order.
    unlikely_indices.sort()

    print("-" * 30)
    print("The indices of the tribes unlikely to be collected are:")

    # Create the final comma-separated string for printing.
    result_string = ", ".join(map(str, unlikely_indices))
    print(result_string)


if __name__ == "__main__":
    analyze_collection_method()