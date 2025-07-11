import sys
# A simple way to suppress warnings in a very controlled manner, if needed.
# For this script, it's not strictly necessary but good practice.
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

def find_latest_battle():
    """
    Finds and prints the latest historical battle depicted in a known painting by Lady Butler.
    """
    # Step 1: Create a dictionary of Lady Butler's paintings and the years of the battles they depict.
    # The format is { "Painting Title": ("Battle/Event Name", Year) }
    paintings = {
        "The Roll Call": ("Crimean War (after the Battle of Inkerman)", 1854),
        "Quatre Bras": ("Battle of Quatre Bras", 1815),
        "The Defence of Rorke's Drift": ("Battle of Rorke's Drift", 1879),
        "Scotland for Ever!": ("Battle of Waterloo", 1815),
        "The Camel Corps": ("Nile Expedition (Battle of Abu Klea)", 1885),
        "The Charge of the Dorset Yeomanry at Agagia": ("Action at Agagia", 1916),
        "Within Sound of the Guns": ("Second Boer War", 1902), # Approximate end year of the war
        "A Detachment of the D.C.L.I. on the March, South Africa": ("Second Boer War", 1900)
    }

    # Step 2: Find the entry with the latest year.
    # We can use the max() function with a lambda key to find the item with the highest year.
    try:
        latest_painting_title, (latest_battle_name, latest_year) = max(paintings.items(), key=lambda item: item[1][1])
    except ValueError:
        print("The list of paintings is empty.")
        return

    # Step 3: Print the result clearly.
    print("Comparing the dates of battles depicted in paintings by Lady Butler...")
    print("-" * 50)
    print(f"The latest battle found is the '{latest_battle_name}'.")
    print(f"The event occurred in the year {latest_year}.")
    print(f"It was depicted in the painting titled: '{latest_painting_title}'.")


# Execute the function
find_latest_battle()