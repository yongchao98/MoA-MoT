def solve_r_flash_champions():
    """
    This function identifies and prints the names of League of Legends champions from a predefined list
    who can perform the R-Flash combo up until Season 14 (2024).
    """
    # The list of champions who can perform the R-Flash combo based on game mechanics.
    # The R ability must have a cast time during which Flash can be used to reposition
    # the ability's origin or area of effect.
    r_flash_champions = [
        "Amumu",
        "Azir",
        "Braum",
        "Cassiopeia",
        "Gnar",
        "Graves",
        "Lee Sin",
        "Neeko",
        "Nilah",
        "Poppy",
        "Qiyana",
        "Rell",
        "Riven",
        "Seraphine",
        "Sona",
        "Thresh",
        "Xin Zhao",
        "Yone"
    ]

    # Joining the list of names into a single string, separated by commas.
    result = ", ".join(r_flash_champions)
    
    # Printing the final result.
    print(result)

solve_r_flash_champions()