def solve():
    """
    This function identifies and prints the names of League of Legends champions
    from a predefined list who can perform the R-Flash combo.
    """
    # List of champions who can perform the R-Flash combo as described.
    # The technique requires a cast time on the ultimate ability that can be
    # buffered with Flash to reposition the ability's effect.
    champions_with_r_flash = [
        "Amumu",
        "Azir",
        "Braum",
        "Cassiopeia",
        "Diana",
        "Evelynn",
        "Gnar",
        "Graves",
        "Illaoi",
        "Jarvan",
        "Lee Sin",
        "Neeko",
        "Nilah",
        "Poppy",
        "Qiyana",
        "Rell",
        "Renata",
        "Riven",
        "Sejuani",
        "Seraphine",
        "Sona",
        "Thresh",
        "Vi",
        "Xin Zhao",
        "Yone"
    ]

    # Print the final list, with names separated by a comma.
    print(",".join(champions_with_r_flash))

solve()