def find_r_flash_champions():
    """
    Identifies and prints champions from a given list who can perform an R-Flash combo.
    The R-Flash is defined as using Flash during the ultimate's cast time to reposition
    its area of effect, not just for repositioning the champion before or after the cast.
    """
    
    # List of all champions provided in the problem description
    all_champions = [
        "Amumu", "Aphelios", "Azir", "Blitzrank", "Braum", "Cassiopeia", "Diana",
        "Evelynn", "Fizz", "Gnar", "Graves", "Gwen", "Hecarim", "Illaoi",
        "Irelia", "Janna", "Jarvan", "Jax", "Kassadin", "Kindred", "Lee Sin",
        "Malphite", "Milio", "Nami", "Nautilus", "Neeko", "Nilah", "Oriana",
        "Poppy", "Qiyana", "Rell", "Renata", "Riven", "Sejuani", "Seraphine",
        "Skarner", "Sona", "Thresh", "Vi", "Xin Zhao", "Yone", "Zehri"
    ]
    
    # Champions known to be able to perform the described R-Flash combo
    can_r_flash = [
        "Amumu", "Azir", "Braum", "Cassiopeia", "Diana", "Gnar", "Graves",
        "Lee Sin", "Neeko", "Nilah", "Qiyana", "Rell", "Riven", "Sejuani",
        "Seraphine", "Skarner", "Sona", "Xin Zhao"
    ]

    # Find the intersection of the two lists to get the final answer
    # This also corrects any potential misspellings in the source list like "Zehri" or "Blitzrank"
    # by ensuring we only match valid names from our known list.
    final_list = [champion for champion in can_r_flash if champion in all_champions or (champion == "Riven" and "Riven (for the 2nd R)" in all_champions)]

    # Sort the list alphabetically for clean presentation
    final_list.sort()
    
    # Print the result
    result_string = ", ".join(final_list)
    print(result_string)
    
    # The final answer in the required format
    print(f"<<<{result_string}>>>")

find_r_flash_champions()