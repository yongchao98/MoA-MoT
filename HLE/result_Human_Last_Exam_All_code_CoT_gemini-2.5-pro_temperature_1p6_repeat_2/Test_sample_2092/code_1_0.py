def solve_r_flash_champions():
    """
    This function identifies League of Legends champions from a predefined list
    who can perform the R-Flash combo up to Season 14 (2024).
    The R-Flash combo is defined as casting the ultimate ability (R) and then
    using Flash during the cast/channel time to buffer the ability,
    changing its origin or effect for optimal damage or crowd control.
    """

    # List of champions that can perform the R-Flash buffer combo
    # based on game mechanics up to and including Season 14 (2024).
    r_flash_champs = [
        "Amumu",
        "Braum",
        "Cassiopeia",
        "Evelynn",
        "Gnar",
        "Graves",
        "Gwen",
        "Irelia",
        "Jax",
        "Lee Sin",
        "Nami",
        "Neeko",
        "Nilah",
        "Orianna",
        "Qiyana",
        "Rell",
        "Riven",
        "Sejuani",
        "Seraphine",
        "Skarner",
        "Sona",
        "Thresh",
        "Xin Zhao",
        "Yone"
    ]

    # Print the final list of champion names, separated by a comma.
    print(",".join(r_flash_champs))

solve_r_flash_champions()