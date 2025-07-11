import pandas as pd
from io import StringIO

def get_latest_battle():
    """
    This function finds the latest historical battle depicted in a painting by Lady Butler.
    """
    # Data of Lady Butler's paintings and the year of the battle depicted.
    data = """Painting,Battle,Year
The Roll Call,Battle of Inkerman,1854
Quatre Bras,Battle of Quatre Bras,1815
Balaclava,Battle of Balaclava,1854
The Defence of Rorke's Drift,Battle of Rorke's Drift,1879
Scotland for Ever!,Battle of Waterloo,1815
Tel-el-Kebir,Battle of Tel El Kebir,1882
To the Front,Franco-Prussian War,1870
Steady the Drums and Fifes,Battle of Albuera,1811
"""

    # Read the data into a pandas DataFrame
    df = pd.read_csv(StringIO(data))

    # Find the painting with the latest battle year
    latest_battle_info = df.loc[df['Year'].idxmax()]

    painting = latest_battle_info['Painting']
    battle = latest_battle_info['Battle']
    year = latest_battle_info['Year']

    print(f"List of battles depicted by Lady Butler and their years:")
    for index, row in df.sort_values('Year').iterrows():
        print(f"- {row['Battle']}: {row['Year']}")
    
    print("\nDetermining the latest battle:")
    print(f"The maximum year found is {year}.")
    print(f"\nThe latest historical battle depicted in a painting by Lady Butler is the '{battle}'.")
    print(f"This event, which took place in the year {year}, was portrayed in her painting titled '{painting}'.")

get_latest_battle()
<<<The Battle of Tel El Kebir>>>