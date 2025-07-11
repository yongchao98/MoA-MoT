import pandas as pd

def find_longest_rivers():
    """
    This function identifies the three longest rivers mentioned in Virgil's Aeneid.
    It uses a predefined dictionary of rivers and their lengths, sorts them,
    and prints the top three.
    """
    # A dictionary of rivers mentioned in the Aeneid and their modern lengths in kilometers.
    # The list includes major rivers whose lengths are significant.
    rivers_in_aeneid = {
        "Nile": 6650,
        "Danube (Ister)": 2850,
        "Euphrates": 2800,
        "Ganges": 2704,
        "Rhine (Rhenus)": 1233,
        "Eridanus (Po)": 652,
        "Tiber": 406,
        "Hebrus (Maritsa)": 480
    }

    # Sort the dictionary by river length in descending order
    sorted_rivers = sorted(rivers_in_aeneid.items(), key=lambda item: item[1], reverse=True)

    # Get the top three longest rivers
    top_three_rivers = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for i, (river, length) in enumerate(top_three_rivers):
        print(f"{i+1}. {river}: approximately {length} km")

if __name__ == "__main__":
    find_longest_rivers()