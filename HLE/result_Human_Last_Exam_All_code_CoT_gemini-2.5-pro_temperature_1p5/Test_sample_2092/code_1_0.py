def find_r_flash_champions():
    """
    This function identifies champions from a predefined list who can perform an R-Flash combo.
    The R-Flash mechanic involves casting the ultimate ability (R) and then using Flash
    during the cast time or channel to reposition the ability's effect, making it more
    effective for damage or crowd control.
    """
    
    champion_list = [
        "Amumu", "Aphelios", "Azir", "Blitzrank", "Braum", "Cassiopeia", "Diana",
        "Evelynn", "Fizz", "Gnar", "Graves", "Gwen", "Hecarim", "Illaoi", "Irelia",
        "Janna", "Jarvan", "Jax", "Kassadin", "Kindred", "Lee Sin", "Malphite",
        "Milio", "Nami", "Nautilus", "Neeko", "Nilah", "Oriana", "Poppy", "Qiyana",
        "Rell", "Renata", "Riven", "Sejuani", "Seraphine", "Skarner", "Sona",
        "Thresh", "Vi", "Xin Zhao", "Yone", "Zehri"
    ]
    
    # Champions who can perform a true, buffered R-Flash for combat effectiveness
    # based on game mechanics up to Season 14 (2024).
    r_flash_champs = [
        "Amumu",
        "Braum",
        "Cassiopeia",
        "Gnar", # As Mega Gnar
        "Lee Sin",
        "Milio",
        "Neeko",
        "Nilah",
        "Poppy", # For uncharged (tap) R
        "Qiyana",
        "Rell",
        "Riven", # For the second R activation (Wind Slash)
        "Sejuani",
        "Seraphine",
        "Skarner",
        "Thresh",
        "Vi",
        "Xin Zhao",
        "Yone"
    ]

    print(",".join(r_flash_champs))

find_r_flash_champions()