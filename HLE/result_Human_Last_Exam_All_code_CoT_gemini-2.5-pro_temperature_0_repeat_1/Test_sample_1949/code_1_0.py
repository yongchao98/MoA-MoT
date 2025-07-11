def find_longest_rivers_in_aeneid():
    """
    This function identifies the three longest rivers mentioned in the Aeneid.
    It stores the rivers and their lengths in a dictionary, sorts them,
    and prints the top three.
    """
    # A dictionary of rivers explicitly mentioned in the Aeneid and their approximate lengths in km.
    # Sources for mentions: Nile (VIII.711), Euphrates (VIII.726), Rhine (VIII.727),
    # Araxes (VIII.728), Po (as Eridanus, VI.659), Hebrus (I.317), Tiber (passim).
    rivers = {
        "Nile": 6650,
        "Euphrates": 2800,
        "Rhine": 1233,
        "Araxes": 1072,
        "Po": 652,
        "Hebrus": 480,
        "Tiber": 406,
    }

    # Sort the rivers by length in descending order
    sorted_rivers = sorted(rivers.items(), key=lambda item: item[1], reverse=True)

    # Get the top three longest rivers
    top_three_rivers = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for river, length in top_three_rivers:
        print(f"- {river}: approximately {length} km")

if __name__ == "__main__":
    find_longest_rivers_in_aeneid()
