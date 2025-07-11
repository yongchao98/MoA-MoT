import pandas as pd

def find_longest_rivers_in_aeneid():
    """
    Identifies the three longest rivers mentioned in the Aeneid.
    
    The function contains a predefined dictionary of rivers mentioned in the text
    and their approximate real-world lengths in kilometers. It then sorts this
    data to find and print the top three longest rivers.
    """
    
    # Dictionary of rivers mentioned in the Aeneid and their approx. lengths (km)
    # Sources include descriptions of places, Aeneas's shield (Book VIII),
    # and travels.
    rivers_in_aeneid = {
        'Nile': 6650,
        'Danube': 2850,
        'Euphrates': 2800,
        'Rhine': 1233,
        'Po': 652,
        'Hebrus': 480,
        'Tiber': 406,
        'Xanthus': 100,
    }

    # Sort the rivers by length in descending order
    # sorted() returns a list of tuples: [('Nile', 6650), ('Danube', 2850), ...]
    sorted_rivers = sorted(rivers_in_aeneid.items(), key=lambda item: item[1], reverse=True)
    
    # Get the top three
    top_three_rivers = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for i, (river, length) in enumerate(top_three_rivers, 1):
        print(f"{i}. {river}: approximately {length} km")

# Execute the function
find_longest_rivers_in_aeneid()