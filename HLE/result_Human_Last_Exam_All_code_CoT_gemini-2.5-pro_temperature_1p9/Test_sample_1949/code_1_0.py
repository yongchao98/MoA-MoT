def find_longest_rivers():
    """
    This function identifies the three longest rivers mentioned in Virgil's Aeneid.
    It stores a list of rivers from the epic and their approximate modern lengths,
    sorts them, and prints the top three.
    """
    # A dictionary of rivers mentioned in the Aeneid and their lengths in km.
    # The list includes major world rivers referenced in the text.
    rivers_in_aeneid = {
        "Nile": 6650,
        "Indus (represented by its tributary Hydaspes)": 3180,
        "Danube (Ister)": 2850,
        "Euphrates": 2800,
        "Ganges": 2525,
        "Tigris": 1850,
        "Rhine": 1233,
        "Araxes": 1072,
        "Po (Eridanus)": 652,
        "Tiber": 406,
    }

    # Sort the dictionary items by length (the value) in descending order.
    sorted_rivers = sorted(rivers_in_aeneid.items(), key=lambda item: item[1], reverse=True)

    # Get the top three longest rivers.
    top_three_rivers = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for i, (river, length) in enumerate(top_three_rivers):
        print(f"{i+1}. {river}: approximately {length} km")

if __name__ == "__main__":
    find_longest_rivers()