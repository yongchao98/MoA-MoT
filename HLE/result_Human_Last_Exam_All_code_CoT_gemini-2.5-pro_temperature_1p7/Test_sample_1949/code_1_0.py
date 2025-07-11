import operator

def find_longest_rivers_in_aeneid():
    """
    Identifies the three longest rivers mentioned in Virgil's Aeneid.

    This function contains a predefined dictionary of rivers mentioned in the text
    and their corresponding lengths. It sorts them to find the top three and
    prints the result.
    """
    # A dictionary of rivers mentioned in the Aeneid and their approximate lengths in kilometers.
    # Latin names are included for reference.
    rivers_in_aeneid = {
        'Nile (Nilus)': 6650,
        'Danube (Ister)': 2850,
        'Euphrates': 2800,
        'Ganges': 2525,
        'Don (Tanais)': 1870,
        'Rhine (Rhenus)': 1233,
        'Po (Padus/Eridanus)': 652,
        'Tiber (Thybris)': 406,
        'Xanthus (Scamander)': 100,
    }

    # Sort the rivers by length in descending order
    sorted_rivers = sorted(rivers_in_aeneid.items(), key=operator.itemgetter(1), reverse=True)

    # Get the top three longest rivers
    top_three_rivers = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for i, (river, length) in enumerate(top_three_rivers):
        print(f"{i+1}. {river}: {length} km")

if __name__ == "__main__":
    find_longest_rivers_in_aeneid()