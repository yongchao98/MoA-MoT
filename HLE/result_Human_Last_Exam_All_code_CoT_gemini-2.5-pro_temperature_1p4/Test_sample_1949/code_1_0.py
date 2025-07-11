import pandas as pd

def find_longest_rivers():
    """
    Identifies the three longest rivers mentioned in the Aeneid.
    """
    # Data of rivers mentioned in the Aeneid and their approximate lengths in km.
    # The Aeneid mentions the Danube by its ancient name, Ister.
    river_data = {
        'Nile': 6650,
        'Danube (Ister)': 2850,
        'Euphrates': 2800,
        'Ganges': 2525,
        'Don (Tanais)': 1870,
        'Tigris': 1850,
        'Rhine': 1233,
        'Po (Eridanus)': 652,
        'Tiber': 406,
        'Simois': 93, # Small river near Troy
        'Scamander (Xanthus)': 100, # Small river near Troy
    }

    # Create a pandas Series for easy sorting
    rivers = pd.Series(river_data)

    # Sort the rivers by length in descending order
    longest_rivers = rivers.sort_values(ascending=False)

    # Get the top three
    top_three = longest_rivers.head(3)

    print("The three longest rivers mentioned in the Aeneid are:")
    for river, length in top_three.items():
        print(f"- {river}: approximately {length} km")

if __name__ == "__main__":
    find_longest_rivers()