import collections

def find_longest_rivers_in_aeneid():
    """
    Identifies the three longest rivers mentioned in the Aeneid
    by their approximate real-world lengths.
    """
    # A dictionary of rivers mentioned in the Aeneid and their approximate lengths in kilometers.
    # The rivers are mentioned in various contexts (geographical, prophetic, mythological).
    # For example:
    # - Nile: Book VIII, on the shield of Aeneas depicting the Battle of Actium.
    # - Danube (Ister): Book VIII, also on the shield.
    # - Euphrates: Book VIII, shown as a conquered river.
    # - Po (Eridanus): Book VI, a river in the underworld.
    # - Tiber: Central to the story's setting in Italy.
    rivers_data = {
        "Nile": 6650,
        "Danube": 2850,
        "Euphrates": 2800,
        "Po": 652,
        "Tiber": 406
        # Other smaller rivers like the Simois, Scamander, and Hebrus are also mentioned
        # but are significantly shorter than the ones listed above.
    }

    # Sort the dictionary items by length (the second element of the tuple) in descending order.
    sorted_rivers = sorted(rivers_data.items(), key=lambda item: item[1], reverse=True)

    # Get the top three.
    top_three = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")

    # Print the result, showing each river's name and its length.
    # The numbers in the final output are the approximate lengths in kilometers.
    first_river_name, first_river_length = top_three[0]
    second_river_name, second_river_length = top_three[1]
    third_river_name, third_river_length = top_three[2]
    
    print(f"1. {first_river_name} (Length: {first_river_length} km)")
    print(f"2. {second_river_name} (Length: {second_river_length} km)")
    print(f"3. {third_river_name} (Length: {third_river_length} km)")

find_longest_rivers_in_aeneid()