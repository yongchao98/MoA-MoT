def solve():
    """
    This function identifies and prints the names of League of Legends champions
    from a predefined list who can perform the R-Flash combo.
    """
    
    # The list of all champions provided in the prompt.
    champion_pool = [
        "Amumu", "Aphelios", "Azir", "Blitzrank", "Braum", "Cassiopeia", 
        "Diana", "Evelynn", "Fizz", "Gnar", "Graves", "Gwen", "Hecarim", 
        "Illaoi", "Irelia", "Janna", "Jarvan", "Jax", "Kassadin", "Kindred", 
        "Lee Sin", "Malphite", "Milio", "Nami", "Nautilus", "Neeko", "Nilah", 
        "Oriana", "Poppy", "Qiyana", "Rell", "Renata", "Riven", "Sejuani", 
        "Seraphine", "Skarner", "Sona", "Thresh", "Vi", "Xin Zhao", "Yone", "Zehri"
    ]

    # The list of champions confirmed to be able to perform the R-Flash combo.
    # Note: "Renata" is used as per the prompt's spelling.
    # Note: "Skarner" is included based on his pre-rework kit, relevant for the "up until season 14" timeframe.
    r_flash_champions = [
        "Amumu",
        "Braum",
        "Cassiopeia",
        "Diana",
        "Gnar",
        "Lee Sin",
        "Nami",
        "Neeko",
        "Nilah",
        "Qiyana",
        "Rell",
        "Renata",
        "Riven",
        "Sejuani",
        "Seraphine",
        "Skarner",
        "Sona",
        "Xin Zhao",
        "Yone"
    ]
    
    # Format the list into a single string, with names separated by ", "
    result = ", ".join(r_flash_champions)
    
    # Print the final formatted string
    print(result)

solve()
<<<Amumu, Braum, Cassiopeia, Diana, Gnar, Lee Sin, Nami, Neeko, Nilah, Qiyana, Rell, Renata, Riven, Sejuani, Seraphine, Skarner, Sona, Xin Zhao, Yone>>>