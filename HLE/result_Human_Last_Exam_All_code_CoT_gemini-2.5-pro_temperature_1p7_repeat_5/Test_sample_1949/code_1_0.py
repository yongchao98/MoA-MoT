def find_longest_aeneid_rivers():
    """
    Identifies and prints the three longest rivers mentioned in the Aeneid.
    
    The script stores a list of rivers mentioned in the epic poem along with their
    approximate modern lengths. It then sorts them to find the three longest.
    The instruction 'output each number in the final equation!' is interpreted
    as clearly displaying the lengths of these rivers.
    """
    
    # Data: A list of dictionaries, each containing a river's name and length in km.
    # These are some of the prominent real-world rivers mentioned in Virgil's Aeneid.
    rivers_in_aeneid = [
        {'name': 'Nile', 'length_km': 6650},
        {'name': 'Danube', 'length_km': 2850},
        {'name': 'Euphrates', 'length_km': 2800},
        {'name': 'Ganges', 'length_km': 2525},
        {'name': 'Tigris', 'length_km': 1850},
        {'name': 'Rhine', 'length_km': 1233},
        {'name': 'Po (Eridanus)', 'length_km': 652},
        {'name': 'Tiber', 'length_km': 406}
    ]

    # Sort the list of rivers by length in descending order
    sorted_rivers = sorted(rivers_in_aeneid, key=lambda x: x['length_km'], reverse=True)

    # Get the top three longest rivers
    top_three = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    
    # Print the name and length of each of the top three rivers
    for river in top_three:
        name = river['name']
        length = river['length_km']
        print(f"{name}: {length} km")

# Execute the function
find_longest_aeneid_rivers()