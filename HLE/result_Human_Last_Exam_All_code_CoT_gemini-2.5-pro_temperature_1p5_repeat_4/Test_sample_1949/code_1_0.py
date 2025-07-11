def find_longest_rivers_in_aeneid():
    """
    This function identifies the three longest rivers mentioned in the Aeneid
    by sorting a predefined list of rivers and their lengths.
    """
    # A dictionary of rivers mentioned in the Aeneid and their approximate lengths in kilometers.
    # Sources for mentions include: Nile (Aen. 8), Danube/Ister (Georgics), Rhine (Aen. 8),
    # Ganges (Aen. 8), Euphrates (Aen. 8), Tiber (frequently), Po/Eridanus (Aen. 6).
    rivers = {
        "Nile": 6650,
        "Danube": 2850,
        "Euphrates": 2800,
        "Ganges": 2525,
        "Tigris": 1900,
        "Rhine": 1233,
        "Po (Eridanus)": 652,
        "Tiber": 406
    }

    # Sort the rivers by length in descending order
    # sorted() returns a list of tuples (key, value)
    sorted_rivers = sorted(rivers.items(), key=lambda item: item[1], reverse=True)

    # Get the top three longest rivers
    top_three = sorted_rivers[:3]

    # Print the result
    print("The three longest rivers mentioned in the Aeneid are:")
    for i, (river, length) in enumerate(top_three):
        print(f"{i+1}. {river}: approximately {length} km")

if __name__ == "__main__":
    find_longest_rivers_in_aeneid()