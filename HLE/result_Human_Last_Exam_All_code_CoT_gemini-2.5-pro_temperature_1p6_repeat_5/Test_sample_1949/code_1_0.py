import pandas as pd

def find_longest_rivers():
    """
    Identifies the three longest rivers mentioned in the Aeneid.
    """
    # Data of rivers mentioned in the Aeneid and their approximate lengths in km.
    # The Aeneid mentions rivers like the Nile (Nilus), Danube (Ister),
    # Euphrates, Tigris, Ganges, Rhine (Rhenus), Po (Eridanus), and Tiber (Tiberis).
    river_data = {
        'Nile': 6650,
        'Danube': 2850,
        'Euphrates': 2800,
        'Ganges': 2525,
        'Tigris': 1850,
        'Rhine': 1230,
        'Po': 652,
        'Tiber': 406
    }

    # Convert the dictionary to a pandas DataFrame for easy sorting.
    df = pd.DataFrame(list(river_data.items()), columns=['River', 'Length (km)'])

    # Sort the DataFrame by length in descending order.
    df_sorted = df.sort_values(by='Length (km)', ascending=False)

    # Get the top three longest rivers.
    top_three_rivers = df_sorted.head(3)

    print("The three longest rivers mentioned in the Aeneid are:")
    for index, row in top_three_rivers.iterrows():
        print(f"- {row['River']}: approximately {row['Length (km)']} km")

find_longest_rivers()