def find_longest_rivers_in_aeneid():
    """
    This function identifies the three longest rivers mentioned in the Aeneid.
    It stores a dictionary of rivers mentioned in the epic and their approximate lengths,
    sorts them, and prints the top three.
    """
    # A dictionary of prominent real-world rivers mentioned in the Aeneid and their lengths in km.
    # Virgil mentions these rivers to illustrate the extent of future Roman power.
    # For example, in Book VIII, Vulcan's shield depicts scenes including the Nile,
    # and in Book VI, the river Eridanus (the Po) is mentioned in the underworld.
    rivers = {
        "Nile": 6650,      # Mentioned in Book VIII
        "Danube": 2850,    # Mentioned as 'Ister' in Book VIII
        "Euphrates": 2800, # Mentioned in Book VIII
        "Ganges": 2525,    # Mentioned in Book IX
        "Rhine": 1233,     # Mentioned in Book VIII
        "Po": 652,         # Mentioned as 'Eridanus' in Book VI
        "Tiber": 406,      # A central river in the latter half of the Aeneid
        "Hebrus": 480      # A river in Thrace, mentioned in Book I
    }

    # Sort the rivers by length in descending order
    # The `sorted` function takes the dictionary's items and sorts them
    # The `key` lambda function tells sorted to use the second element of each item (the length)
    # `reverse=True` ensures the sort is from longest to shortest
    sorted_rivers = sorted(rivers.items(), key=lambda item: item[1], reverse=True)

    # Get the top three longest rivers
    top_three = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    
    # Print each of the top three rivers with their name and length
    for i, (river, length) in enumerate(top_three, 1):
        print(f"{i}. {river}: {length} km")

# Execute the function
find_longest_rivers_in_aeneid()