def find_longest_rivers():
    """
    This function identifies the three longest rivers mentioned in Virgil's Aeneid.
    The data is pre-compiled from literary and geographical sources.
    """
    
    # Data: A list of prominent rivers mentioned in the Aeneid and their approximate lengths in km.
    # The rivers and their mentions are:
    # Nile: On Aeneas's shield (Book VIII)
    # Danube: Mentioned as 'Ister' (Book III)
    # Euphrates: On Aeneas's shield (Book VIII)
    # Tigris: (Book VIII)
    # Po: Mentioned as 'Eridanus' (Book VI)
    # Hebrus: Thracian river (Book I)
    # Tiber: The river of Rome, central to the story (multiple books)
    # Scamander: River of Troy (multiple books)
    
    rivers_in_aeneid = [
        ("Nile", 6650),
        ("Danube", 2850),
        ("Euphrates", 2800),
        ("Tigris", 1900),
        ("Po", 652),
        ("Hebrus", 480),
        ("Tiber", 406),
        ("Scamander", 110)
    ]

    # Sort the rivers by length in descending order.
    # The `key=lambda item: item[1]` tells the sort function to use the second element (length) for comparison.
    sorted_rivers = sorted(rivers_in_aeneid, key=lambda item: item[1], reverse=True)
    
    # Get the top three longest rivers.
    top_three_rivers = sorted_rivers[:3]
    
    print("The three longest rivers mentioned in the Aeneid are:")
    
    # Print the name and length of each of the top three rivers.
    for i, river in enumerate(top_three_rivers):
        # The equation here is simply Rank = River: Length
        print(f"{i + 1}. {river[0]}: {river[1]} km")

if __name__ == "__main__":
    find_longest_rivers()