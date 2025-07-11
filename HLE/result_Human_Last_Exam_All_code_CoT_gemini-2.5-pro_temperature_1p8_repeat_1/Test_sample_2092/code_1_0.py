def solve_r_flash_puzzle():
    """
    Identifies and prints the names of League of Legends champions from a predefined list
    who can perform an R-Flash combo to buffer their ultimate ability.
    """

    # The list of champions from the problem description who can perform a buffered R-Flash.
    # This is based on in-game mechanics where pressing R then Flash repositions the ability's effect.
    r_flash_champions = [
        "Amumu",
        "Azir",
        "Braum",
        "Cassiopeia",
        "Evelynn",
        "Gnar",
        "Graves",
        "Illaoi",
        "Janna",
        "Lee Sin",
        "Nami",
        "Neeko",
        "Nilah",
        "Oriana",
        "Poppy",
        "Qiyana",
        "Rell",
        "Renata",
        "Sejuani",
        "Seraphine",
        "Skarner",
        "Sona",
        "Thresh",
        "Xin Zhao",
        "Yone"
    ]
    
    # The prompt asks for Riven's 2nd R, which cannot be buffered this way.
    # Therefore, Riven is excluded from the final list.
    
    # We join the names with ", " to create the final output string.
    result_string = ", ".join(r_flash_champions)
    
    print(result_string)

solve_r_flash_puzzle()