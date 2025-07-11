import operator

def find_latest_battle():
    """
    This function identifies the latest historical battle depicted in a painting by Lady Butler.
    """
    # A dictionary of battles/campaigns depicted in Lady Butler's paintings and their years.
    battles = {
        "Battle of Albuera ('Steady the Drums and Fifes')": 1811,
        "Battle of Waterloo ('Scotland for Ever!')": 1815,
        "Battle of Quatre Bras ('Quatre Bras')": 1815,
        "Battle of Inkerman ('The Roll Call')": 1854,
        "Battle of Rorke's Drift ('The Defence of Rorke's Drift')": 1879,
        "First Boer War / Laing's Nek ('Floreat Etona!')": 1881,
        "Battle of Tel el-Kebir ('After the Battle')": 1882,
        "Gordon Relief Expedition ('A Desert Grave')": 1885,
    }

    print("Comparing the years of events depicted in Lady Butler's paintings:")
    # Print each battle and its year for comparison
    for battle, year in battles.items():
        print(f"{battle}: {year}")

    # Find the battle with the latest year
    latest_battle = max(battles.items(), key=operator.itemgetter(1))
    
    # Print the equation representing the comparison
    # Example: max(1811, 1815, 1815, 1854, 1879, 1881, 1882, 1885) = 1885
    year_list = list(battles.values())
    print(f"\nFinding the maximum year in the list: max({', '.join(map(str, year_list))}) = {latest_battle[1]}")
    
    print(f"\nThe latest historical event depicted in a painting by Lady Butler is the {latest_battle[0].split('(')[0].strip()} which took place in {latest_battle[1]}.")

find_latest_battle()