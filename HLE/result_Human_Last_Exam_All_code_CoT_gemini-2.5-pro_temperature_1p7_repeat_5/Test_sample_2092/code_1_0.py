def solve_r_flash():
    """
    Identifies and prints the names of champions from a predefined list
    who can perform the R-Flash combo in League of Legends.
    """
    # List of all champions provided by the user.
    all_champions = [
        "Amumu", "Aphelios", "Azir", "Blitzrank", "Braum", "Cassiopeia",
        "Diana", "Evelynn", "Fizz", "Gnar", "Graves", "Gwen", "Hecarim",
        "Illaoi", "Irelia", "Janna", "Jarvan", "Jax", "Kassadin", "Kindred",
        "Lee Sin", "Malphite", "Milio", "Nami", "Nautilus", "Neeko", "Nilah",
        "Oriana", "Poppy", "Qiyana", "Rell", "Renata", "Riven", "Sejuani",
        "Seraphine", "Skarner", "Sona", "Thresh", "Vi", "Xin Zhao", "Yone",
        "Zehri"
    ]
    
    # Champions confirmed to have a true R-Flash combo mechanic.
    r_flash_champions = [
        "Amumu", "Azir", "Braum", "Cassiopeia", "Evelynn", "Fizz", "Gnar", 
        "Graves", "Hecarim", "Lee Sin", "Neeko", "Nilah", "Oriana", "Qiyana", 
        "Rell", "Riven", "Sejuani", "Seraphine", "Sona", "Thresh", "Xin Zhao", 
        "Yone"
    ]

    # The user asked for Riven (for the 2nd R), but we'll just list "Riven".
    # The user misspelled "Blitzcrank" as "Blitzrank" and "Zeri" as "Zehri".
    # The code will check against the canonical list of R-Flash champs.
    
    # Filter the original list to ensure we only include champions from that list.
    # This also preserves the user's original spelling if it matches.
    # For this task, we will just construct the final string from the confirmed list.
    
    final_answer_string = ",".join(r_flash_champions)
    
    print(final_answer_string)

solve_r_flash()
<<<Amumu,Azir,Braum,Cassiopeia,Evelynn,Fizz,Gnar,Graves,Hecarim,Lee Sin,Neeko,Nilah,Oriana,Qiyana,Rell,Riven,Sejuani,Seraphine,Sona,Thresh,Xin Zhao,Yone>>>