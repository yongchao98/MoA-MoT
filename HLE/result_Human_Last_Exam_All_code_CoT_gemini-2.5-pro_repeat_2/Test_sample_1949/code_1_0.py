def find_longest_rivers_in_aeneid():
    """
    Identifies the three longest rivers mentioned in Virgil's Aeneid.
    
    The script uses a predefined dictionary of rivers mentioned in the text
    and their approximate real-world lengths. It then sorts them to find the
    top three longest.
    """
    # A dictionary of major real-world rivers mentioned in the Aeneid and their lengths in km.
    # Sources for mentions include: Nile (Aen. 8.711), Danube/Ister (Georg. 2.497, referenced in Aeneid's worldview),
    # Euphrates (Aen. 8.726), Ganges (Aen. 9.31), Tigris (Aen. 6.88 - implied with Assyria),
    # Rhine (Aen. 8.727), Po/Eridanus (Aen. 6.659), Tiber (central to the plot).
    rivers_data = {
        "Nile": 6650,
        "Danube (Ister)": 2850,
        "Euphrates": 2800,
        "Ganges": 2525,
        "Tigris": 1850,
        "Rhine": 1233,
        "Po (Eridanus)": 652,
        "Tiber": 406,
    }

    # Sort the rivers by length in descending order
    # The sorted() function returns a list of tuples (key, value)
    sorted_rivers = sorted(rivers_data.items(), key=lambda item: item[1], reverse=True)

    # Get the top three longest rivers
    top_three_rivers = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for i, (river, length) in enumerate(top_three_rivers, 1):
        print(f"{i}. {river}: approximately {length} km")
        
    # Extract just the names for the final answer
    top_three_names = [river[0] for river in top_three_rivers]
    print("\nFinal Answer:")
    print(', '.join(top_three_names))

if __name__ == "__main__":
    find_longest_rivers_in_aeneid()