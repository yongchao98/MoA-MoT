def solve_r_flash_champions():
    """
    This function identifies and prints the names of League of Legends champions
    from a given list who can perform a true R-Flash combo up to Season 14.
    A true R-Flash is defined as using Flash during the ultimate's cast time
    to buffer the ability and alter its effective area for offensive or
    crowd-control purposes.
    """

    champions_who_can_r_flash = [
        "Amumu",
        "Braum",
        "Cassiopeia",
        "Diana",
        "Gnar",
        "Hecarim",
        "Lee Sin",
        "Neeko",
        "Nilah",
        "Qiyana",
        "Rell",
        "Riven", # As specified, this applies to her 2nd R cast (Wind Slash)
        "Sejuani",
        "Seraphine",
        "Skarner",
        "Sona",
        "Xin Zhao",
        "Yone",
        "Zeri" # Name corrected from "Zehri"
    ]

    # Join the list of champion names into a single string, separated by commas
    result = ", ".join(champions_who_can_r_flash)
    
    # Print the final result
    print(result)

solve_r_flash_champions()