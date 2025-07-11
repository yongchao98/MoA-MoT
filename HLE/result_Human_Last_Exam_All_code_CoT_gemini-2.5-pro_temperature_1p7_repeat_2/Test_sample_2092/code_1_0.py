def solve_lol_r_flash():
    """
    Identifies and prints the names of League of Legends champions from a given list
    who could perform an R-Flash combo up until Season 14 (2024).

    The R-Flash combo is defined as buffering Flash during the cast time of the ultimate
    ability (R) to alter its point of origin or area of effect.
    """
    
    # Full list of champions provided by the user.
    all_champions = [
        "Amumu", "Aphelios", "Azir", "Blitzcrank", "Braum", "Cassiopeia", 
        "Diana", "Evelynn", "Fizz", "Gnar", "Graves", "Gwen", "Hecarim", 
        "Illaoi", "Irelia", "Janna", "Jarvan IV", "Jax", "Kassadin", 
        "Kindred", "Lee Sin", "Malphite", "Milio", "Nami", "Nautilus", 
        "Neeko", "Nilah", "Orianna", "Poppy", "Qiyana", "Rell", "Renata", 
        "Riven", "Sejuani", "Seraphine", "Skarner", "Sona", "Thresh", 
        "Vi", "Xin Zhao", "Yone", "Zeri" # Corrected from Zehri
    ]

    # Champions who can perform the R-Flash combo based on game mechanics.
    r_flash_champions = [
        "Amumu",
        "Azir",
        "Braum",
        "Cassiopeia",
        "Gnar",
        "Hecarim",
        "Janna",
        "Jarvan IV", # Name from prompt is "Jarvan"
        "Lee Sin",
        "Neeko",
        "Poppy",
        "Qiyana",
        "Rell",
        "Riven", # As specified, for the 2nd R
        "Sejuani",
        "Seraphine",
        "Skarner",
        "Sona",
        "Thresh",
        "Vi",
        "Xin Zhao",
        "Yone"
    ]
    
    # Correcting "Jarvan" from the prompt to the full name for clarity,
    # but will stick to the prompt's provided names where possible.
    # The name "Jarvan" is used in the prompt instead of "Jarvan IV". Let's use the one from the prompt.
    r_flash_champions[7] = "Jarvan"

    print(", ".join(r_flash_champions))

solve_lol_r_flash()