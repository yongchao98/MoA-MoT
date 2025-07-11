def find_longest_aeneid_rivers():
    """
    Identifies the three longest rivers mentioned in Virgil's Aeneid,
    sorts them by length, and prints the result.
    """
    # A list of tuples containing rivers mentioned in the Aeneid and their approx. lengths (km).
    # Sources include descriptions of Aeneas's shield and other geographical references.
    rivers = [
        ("Tiber", 406),
        ("Po", 652),
        ("Nile", 6650),
        ("Euphrates", 2800),
        ("Rhine", 1233),
        ("Danube", 2850),
        ("Xanthus", 100),
        ("Ganges", 2525),
        ("Indus", 3180),
        ("Araxes", 1072)
    ]

    # Sort the list of rivers by length in descending order.
    # The lambda function specifies that the sorting key is the second element (index 1) of each tuple.
    sorted_rivers = sorted(rivers, key=lambda item: item[1], reverse=True)

    # Select the top three longest rivers from the sorted list.
    top_three = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    
    # Print each of the top three rivers with its rank, name, and length.
    # This loop outputs each number in the final result.
    for i, river in enumerate(top_three):
        name = river[0]
        length = river[1]
        # The rank is the loop index (i) + 1.
        # The name is the first element of the tuple.
        # The length is the second element of the tuple.
        print(f"Rank {i + 1}: The {name} river is {length} km long.")

if __name__ == '__main__':
    find_longest_aeneid_rivers()