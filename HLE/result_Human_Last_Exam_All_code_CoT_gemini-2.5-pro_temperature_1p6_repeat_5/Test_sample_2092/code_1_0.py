def solve_r_flash_champions():
    """
    This function identifies and prints the names of League of Legends champions
    from a predefined list who can perform an R-Flash combo.
    """
    # The list of champions identified as being able to perform the R-Flash combo
    # based on their ability mechanics up until Season 14 (2024).
    # The combo involves using Flash during the ultimate's cast time to reposition its effect.
    r_flash_champs = [
        "Amumu",
        "Azir",
        "Braum",
        "Cassiopeia",
        "Diana",
        "Evelynn",
        "Gnar",
        "Graves",
        "Jax",
        "Lee Sin",
        "Neeko",
        "Qiyana",
        "Rell",
        "Riven",
        "Sejuani",
        "Seraphine",
        "Skarner",
        "Sona",
        "Xin Zhao",
        "Yone"
    ]

    print(", ".join(r_flash_champs))

solve_r_flash_champions()