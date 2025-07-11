import pandas as pd

def find_longest_rivers_aeneid():
    """
    This function identifies the three longest rivers mentioned in the Aeneid.
    It stores a list of rivers mentioned in the text and their approximate lengths,
    sorts them, and prints the top three.
    """
    # Data of rivers mentioned in the Aeneid and their approximate lengths in km.
    # Sources for mentions: Vergil's Aeneid.
    # Sources for length: Geographical atlases and encyclopedias.
    river_data = {
        'River': ['Nile', 'Danube', 'Euphrates', 'Ganges', 'Tigris', 'Rhine', 'Araxes', 'Po', 'Tiber'],
        'Latin Name': ['Nilus', 'Ister', 'Euphrates', 'Ganges', 'Tigris', 'Rhenus', 'Araxes', 'Eridanus', 'Tiberis'],
        'Approx. Length (km)': [6650, 2850, 2800, 2525, 1900, 1230, 1072, 652, 406]
    }

    # Create a DataFrame and sort it by length in descending order
    df = pd.DataFrame(river_data)
    df_sorted = df.sort_values(by='Approx. Length (km)', ascending=False)

    # Get the top three longest rivers
    top_three_rivers = df_sorted.head(3)

    print("The three longest rivers mentioned in the Aeneid are:")
    
    # Print the results one by one
    for i, row in enumerate(top_three_rivers.itertuples(index=False), 1):
        river_name = row[0]
        length = row[2]
        print(f"{i}. {river_name}: {length} km")

# Execute the function
find_longest_rivers_aeneid()