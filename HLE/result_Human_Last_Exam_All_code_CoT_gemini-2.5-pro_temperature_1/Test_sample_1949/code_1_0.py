def find_longest_aeneid_rivers():
    """
    Identifies and prints the three longest rivers mentioned in the Aeneid.
    """
    # A list of tuples containing major rivers mentioned in the Aeneid and their lengths in km.
    # Format: (River Name, Length in km)
    rivers_in_aeneid = [
        ("Nile", 6650),
        ("Danube", 2850),
        ("Euphrates", 2800),
        ("Tigris", 1850),
        ("Rhine", 1233),
        ("Po", 652),
        ("Tiber", 406)
    ]

    # Sort the list of rivers by length in descending order.
    # The 'key' argument uses a lambda function to specify sorting by the second element of each tuple (the length).
    sorted_rivers = sorted(rivers_in_aeneid, key=lambda x: x[1], reverse=True)

    # Get the top three longest rivers from the sorted list.
    top_three = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for i, (river, length) in enumerate(top_three, 1):
        # The final equation is the rank, river name, and its length.
        print(f"{i}. {river}: {length} km")

if __name__ == '__main__':
    find_longest_aeneid_rivers()