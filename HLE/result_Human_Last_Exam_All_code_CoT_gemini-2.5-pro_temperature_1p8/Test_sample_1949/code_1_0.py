def find_longest_aeneid_rivers():
    """
    Identifies the three longest rivers mentioned in Virgil's Aeneid.

    This function contains a dictionary of major real-world rivers cited in the Aeneid
    and their approximate lengths in kilometers. It then sorts these rivers
    to find and print the top three longest ones.
    """
    # A dictionary mapping river names to their approximate lengths in km.
    # These rivers are confirmed to be mentioned in the Aeneid.
    rivers_data = {
        "Nile": 6650,
        "Danube": 2850,
        "Euphrates": 2800,
        "Ganges": 2704,
        "Rhine": 1233,
        "Rhone": 813,
        "Po": 652,
        "Tiber": 405,
        "Hebrus": 480,
        "Maeander": 548
    }

    # Sort the dictionary items based on length (the value) in descending order.
    # The result is a list of (river_name, length) tuples.
    sorted_rivers = sorted(rivers_data.items(), key=lambda item: item[1], reverse=True)

    # Get the top three rivers from the sorted list.
    top_three = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")

    # Print the name and length for each of the top three rivers.
    # For the output "equation", we will explicitly state each number.
    print(f"1. {top_three[0][0]}: approximately {top_three[0][1]} km long.")
    print(f"2. {top_three[1][0]}: approximately {top_three[1][1]} km long.")
    print(f"3. {top_three[2][0]}: approximately {top_three[2][1]} km long.")

if __name__ == '__main__':
    find_longest_aeneid_rivers()