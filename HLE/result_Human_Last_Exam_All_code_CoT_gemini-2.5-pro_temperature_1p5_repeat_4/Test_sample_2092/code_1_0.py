def solve_r_flash_champions():
    """
    This function identifies and prints the names of League of Legends champions from a given list
    who could perform a buffered R-Flash combo up until Season 14 (2024).
    The R-Flash is defined as casting the ultimate ability (R) and then using Flash
    during the cast animation to change the ability's origin or area of effect.
    """
    
    # List of champions that can perform the R-Flash combo based on game mechanics.
    # The ultimate ability has a cast time or wind-up that can be buffered with Flash.
    r_flash_champions = [
        "Amumu",
        "Azir",
        "Braum",
        "Cassiopeia",
        "Diana",
        "Gnar",
        "Graves",
        "Hecarim",
        "Illaoi",
        "Lee Sin",
        "Nami",
        "Neeko",
        "Nilah",
        "Orianna",
        "Poppy",
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
    
    # Joining the list of names with a comma and a space for clean output.
    result = ", ".join(r_flash_champions)
    
    print(result)

solve_r_flash_champions()