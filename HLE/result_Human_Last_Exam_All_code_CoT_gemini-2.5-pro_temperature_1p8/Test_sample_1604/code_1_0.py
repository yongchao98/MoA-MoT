import datetime

def find_latest_battle_painting():
    """
    Finds the latest historical battle depicted in a known painting by Lady Butler.
    """
    # A dictionary mapping painting titles to their depicted battle and date.
    # Format: { "Painting Title": ("Battle Name", datetime.date(YYYY, M, D)) }
    paintings = {
        "Steady the Drums and Fifes": ("Battle of Albuera", datetime.date(1811, 5, 16)),
        "The 28th Regiment at Quatre Bras": ("Battle of Quatre Bras", datetime.date(1815, 6, 16)),
        "Scotland for Ever!": ("Battle of Waterloo", datetime.date(1815, 6, 18)),
        "The Roll Call": ("Battle of Inkerman", datetime.date(1854, 11, 5)),
        "Balaclava": ("Battle of Balaclava", datetime.date(1854, 10, 25)),
        "The Defence of Rorke's Drift": ("Battle of Rorke's Drift", datetime.date(1879, 1, 22)),
        "Tel-el-Kebir": ("Battle of Tel El Kebir", datetime.date(1882, 9, 13)),
        "The Charge of the Dorset Yeomanry at Agagia, 26th February 1916": ("Affair of Agagia", datetime.date(1916, 2, 26)),
    }

    # Use the max() function to find the entry with the latest date.
    # The `key` argument specifies that we should compare the date part of each dictionary value.
    if not paintings:
        print("No paintings found in the dataset.")
        return

    latest_painting_title, (battle_name, battle_date) = max(paintings.items(), key=lambda item: item[1][1])

    # Print the results
    print(f"The latest battle depicted by Lady Butler is the '{battle_name}'.")
    print(f"It is shown in her painting '{latest_painting_title}'.")
    print(f"The date of this battle is: Year={battle_date.year}, Month={battle_date.month}, Day={battle_date.day}")


if __name__ == '__main__':
    find_latest_battle_painting()