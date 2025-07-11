def find_longest_rivers_in_aeneid():
    """
    This function identifies the three longest rivers mentioned in Virgil's Aeneid.
    It stores a list of these rivers and their approximate lengths, sorts them,
    and prints the top three.
    """
    # A list of tuples containing (River Name, Length in km)
    # Sources for mentions are in Virgil's Aeneid, primarily Books VI and VIII.
    rivers = [
        ("Nile", 6650),
        ("Danube", 2850),
        ("Euphrates", 2800),
        ("Tigris", 1900),
        ("Rhine", 1233),
        ("Po", 652),
        ("Tiber", 406)
    ]

    # Sort the list of rivers by length in descending order
    # The key=lambda item: item[1] tells sorted() to use the second element (length) for sorting
    sorted_rivers = sorted(rivers, key=lambda item: item[1], reverse=True)

    # Get the top three longest rivers
    top_three = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for i, river in enumerate(top_three):
        # The output format is "1. River: Length km"
        print(f"{i + 1}. {river[0]}: {river[1]} km")

if __name__ == "__main__":
    find_longest_rivers_in_aeneid()