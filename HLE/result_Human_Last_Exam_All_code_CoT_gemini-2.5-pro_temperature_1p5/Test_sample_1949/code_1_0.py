def find_longest_rivers_in_aeneid():
    """
    Identifies the three longest rivers mentioned in Virgil's Aeneid.

    The function stores a dictionary of major rivers mentioned in the text
    and their approximate real-world lengths in kilometers. It then sorts
    these rivers by length to find the top three.
    """
    # A dictionary of rivers mentioned in the Aeneid and their approx. lengths (km).
    # Sources for mentions include:
    # Nile (Nilus): Aeneid, Book 8, line 711
    # Danube (Ister): Aeneid, Book 8, line 727
    # Euphrates: Aeneid, Book 8, line 726
    # Rhine (Rhenus): Aeneid, Book 8, line 727
    # Po (Eridanus): Aeneid, Book 6, line 659
    # Tiber (Tiberis): Mentioned extensively in the second half of the epic.
    
    river_data = {
        'Nile': 6650,
        'Danube': 2850,
        'Euphrates': 2800,
        'Rhine': 1230,
        'Po': 652,
        'Tiber': 406
    }

    # Sort the dictionary items by length (the second element of the tuple) in descending order
    sorted_rivers = sorted(river_data.items(), key=lambda item: item[1], reverse=True)

    # Get the top three longest rivers
    top_three = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid, based on their real-world lengths, are:")
    for i, (river, length) in enumerate(top_three):
        print(f"{i+1}. {river} (Length: approximately {length} km)")

find_longest_rivers_in_aeneid()