import sys

# Suppress warnings that are not relevant for this simple script
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

def find_latest_battle():
    """
    Finds the latest historical battle depicted in a painting by Lady Butler from a predefined list.
    """
    # A dictionary mapping battle names (and their paintings) to their year of occurrence.
    # Format: "Painting Title (Battle)" : Year
    battles = {
        "'The Roll Call' (Battle of Inkerman)": 1854,
        "'Quatre Bras' (Battle of Quatre Bras)": 1815,
        "'Balaclava' (Battle of Balaclava)": 1854,
        "'The Defence of Rorke's Drift' (Battle of Rorke's Drift)": 1879,
        "'Scotland for Ever!' (Battle of Waterloo)": 1815,
        "'After the Battle' (Battle of Tel-el-Kebir)": 1882,
        "'The Camel Corps' (Battle of Abu Klea)": 1885,
        "'Floreat Etona!' (Battle of Laing's Nek)": 1881,
        "'A Lament After the Battle' (Battle of Colenso)": 1899,
        "'The charge of the 21st Lancers at Omdurman' (Battle of Omdurman)": 1898,
        "'Steady the Drums and Fifes' (Battle of Albuera)": 1811,
        "'The Last Stand of the 44th Foot at Gundamuck' (Retreat from Kabul)": 1842
    }

    # Find the key (painting/battle) with the maximum value (year)
    latest_battle_key = max(battles, key=battles.get)
    latest_year = battles[latest_battle_key]

    print(f"The latest historical battle depicted in a painting by Lady Butler from the list is:")
    print(f"Painting/Battle: {latest_battle_key}")
    print(f"Year: {latest_year}")

find_latest_battle()