def solve():
    """
    This function identifies and prints the names of League of Legends champions
    from a predefined list who can perform the R-Flash combo.
    """
    champions_who_can_r_flash = [
        "Amumu",
        "Braum",
        "Cassiopeia",
        "Diana",
        "Gnar",
        "Lee Sin",
        "Neeko",
        "Orianna",
        "Qiyana",
        "Rell",
        "Riven (for the 2nd R)",
        "Sejuani",
        "Seraphine",
        "Sona",
        "Xin Zhao",
        "Yone"
    ]

    # Join the list of names into a single string, separated by ","
    result = ", ".join(champions_who_can_r_flash)
    
    print(result)

solve()
<<<Amumu, Braum, Cassiopeia, Diana, Gnar, Lee Sin, Neeko, Orianna, Qiyana, Rell, Riven (for the 2nd R), Sejuani, Seraphine, Sona, Xin Zhao, Yone>>>