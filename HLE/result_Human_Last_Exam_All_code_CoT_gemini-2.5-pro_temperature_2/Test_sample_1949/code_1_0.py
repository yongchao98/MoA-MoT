def find_longest_rivers():
    """
    Identifies the three longest rivers mentioned in the Aeneid and prints them.
    """
    # A dictionary of rivers mentioned in the Aeneid and their approximate modern lengths in kilometers.
    # Sources for mentions include the full text of the Aeneid.
    # For example: Nile (Aen. 6.800, 8.711-13), Danube/Ister (Aen. 6.800), Euphrates (Aen. 8.726),
    # Ganges (Aen. 9.31), Tiber (frequent), Rhine/Rhenus (Aen. 8.727), Tigris (Aen. 8.728),
    # and Eridanus/Po (Aen. 6.659).
    river_lengths = {
        "Nile": 6650,
        "Danube (Ister)": 2850,
        "Euphrates": 2800,
        "Ganges": 2525,
        "Tigris": 1850,
        "Rhine (Rhenus)": 1233,
        "Eridanus (Po)": 652,
        "Hebrus": 480,
        "Tiber": 406
        # Other rivers like Simois, Xanthus, Numicus, etc., are much shorter and omitted for clarity.
    }

    # Sort the dictionary items by length (the value) in descending order.
    sorted_rivers = sorted(river_lengths.items(), key=lambda item: item[1], reverse=True)

    # Get the top three longest rivers.
    top_three_rivers = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for i, (name, length) in enumerate(top_three_rivers):
        print(f"{i + 1}. {name}: approximately {length} km")

find_longest_rivers()