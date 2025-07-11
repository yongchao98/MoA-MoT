def find_r_flash_champions():
    """
    This function identifies and prints the names of League of Legends champions
    from a predefined list who can perform a buffered R-Flash combo.
    The list is based on champion mechanics up to Season 14 (2024).
    """
    
    # List of champions who can perform the specified R-Flash combo
    # The combo is defined as casting R then Flashing during the cast/channel time
    # to redirect the ability's damage or crowd-control effect.
    r_flash_champions = [
        "Amumu",
        "Azir",
        "Braum",
        "Cassiopeia",
        "Diana",
        "Gnar",
        "Janna",
        "Lee Sin",
        "Neeko",
        "Nilah",
        "Qiyana",
        "Riven (for the 2nd R)",
        "Sejuani",
        "Seraphine",
        "Skarner",
        "Sona",
        "Xin Zhao",
        "Yone",
        "Zeri" # Correcting the user's typo from "Zehri"
    ]
    
    # Print the final list of champion names separated by commas
    print(", ".join(r_flash_champions))

find_r_flash_champions()
<<<Amumu, Azir, Braum, Cassiopeia, Diana, Gnar, Janna, Lee Sin, Neeko, Nilah, Qiyana, Riven (for the 2nd R), Sejuani, Seraphine, Skarner, Sona, Xin Zhao, Yone, Zeri>>>