def find_longest_aeneid_rivers():
    """
    This script identifies the three longest rivers mentioned in Virgil's Aeneid.
    It uses a predefined dictionary of rivers mentioned in the text and their approximate modern lengths.
    """
    # A dictionary of major rivers mentioned in the Aeneid and their lengths in km.
    # Mythological rivers are excluded.
    rivers_data = {
        "Nile": 6650,
        "Danube": 2850,
        "Euphrates": 2800,
        "Ganges": 2525,
        "Rhine": 1230,
        "Tagus": 1007,
        "Po": 652,
        "Tiber": 406
    }

    # Sort the rivers by length in descending order.
    # The key for sorting is the length (the second item in each key-value pair).
    sorted_rivers = sorted(rivers_data.items(), key=lambda item: item[1], reverse=True)

    # Get the top three longest rivers.
    top_three = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    # Loop through the top three and print the name and length for each.
    # The "equation" here is showing the name and its corresponding value.
    for river, length in top_three:
        print(f"1. {river}: {length} km")
        # Just to fulfill the "final equation" part, let's break down the next two.
        if top_three.index((river, length)) == 0:
            second_river, second_length = top_three[1]
            third_river, third_length = top_three[2]
            print(f"2. {second_river}: {second_length} km")
            print(f"3. {third_river}: {third_length} km")
            break

find_longest_aeneid_rivers()