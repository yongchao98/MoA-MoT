import operator

def find_latest_battle():
    """
    Identifies the latest historical battle from a list of battles depicted
    in paintings by Lady Butler.
    """
    # A dictionary of battles depicted by Lady Butler and their years.
    # Keys are battle names (and associated paintings), values are the years.
    battles = {
        "'The 28th Regiment at Quatre Bras' (Battle of Quatre Bras)": 1815,
        "'Scotland for Ever!' (Battle of Waterloo)": 1815,
        "'The Remnants of an Army' (Retreat from Kabul)": 1842,
        "'The Roll Call' (Battle of Inkerman, Crimean War)": 1854,
        "'Balaclava' (Battle of Balaclava, Crimean War)": 1854,
        "'The Defence of Rorke's Drift' (Battle of Rorke's Drift)": 1879,
        "'Tel-el-Kebir' (Battle of Tel El Kebir)": 1882,
    }

    # Find the battle with the latest year.
    # The max function with operator.itemgetter(1) finds the key-value pair with the highest value (the year).
    latest_battle, year = max(battles.items(), key=operator.itemgetter(1))

    # Print the result.
    print("List of battles and their years depicted by Lady Butler:")
    for battle, y in battles.items():
        print(f"- {battle}: {y}")
    
    print("\n---")
    print(f"The latest historical battle depicted in a painting by Lady Butler is the {latest_battle}, which took place in the year {year}.")

if __name__ == "__main__":
    find_latest_battle()