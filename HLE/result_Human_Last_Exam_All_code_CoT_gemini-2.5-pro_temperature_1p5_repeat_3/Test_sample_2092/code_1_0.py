def find_r_flash_champions():
    """
    Identifies and prints the names of champions from a predefined list who can perform an R-Flash combo.
    The list of champions is based on game mechanics up to season 14 (2024).
    """
    
    # List of champions who can perform a buffered R-Flash.
    # This combo involves using R and then flashing during the cast time
    # to have the ability's effect originate from the new location.
    r_flash_champs = [
        "Amumu",
        "Azir",
        "Braum",
        "Cassiopeia",
        "Gnar",
        "Graves",
        "Hecarim",
        "Jarvan",
        "Lee Sin",
        "Milio",
        "Neeko",
        "Nilah",
        "Oriana",
        "Poppy",
        "Qiyana",
        "Rell",
        "Riven (for the 2nd R)",
        "Seraphine",
        "Skarner",
        "Sona",
        "Thresh",
        "Xin Zhao",
        "Yone",
        "Zehri" # Zeri
    ]
    
    # Print the final list as a comma-separated string
    print(", ".join(r_flash_champs))

find_r_flash_champions()
<<<Amumu, Azir, Braum, Cassiopeia, Gnar, Graves, Hecarim, Jarvan, Lee Sin, Milio, Neeko, Nilah, Oriana, Poppy, Qiyana, Rell, Riven (for the 2nd R), Seraphine, Skarner, Sona, Thresh, Xin Zhao, Yone, Zehri>>>