def solve():
    """
    This function identifies and prints the names of League of Legends champions
    from a predefined list who can perform an R-Flash combo.
    """
    champions_who_can_r_flash = [
        "Amumu",
        "Azir",
        "Braum",
        "Cassiopeia",
        "Gnar",
        "Graves",
        "Hecarim",
        "Jarvan",
        "Lee Sin",
        "Neeko",
        "Oriana",
        "Qiyana",
        "Rell",
        "Riven (for the 2nd R)",
        "Sejuani",
        "Seraphine",
        "Skarner",
        "Sona",
        "Thresh",
        "Xin Zhao",
        "Yone"
    ]
    
    result = ", ".join(champions_who_can_r_flash)
    print(result)

solve()