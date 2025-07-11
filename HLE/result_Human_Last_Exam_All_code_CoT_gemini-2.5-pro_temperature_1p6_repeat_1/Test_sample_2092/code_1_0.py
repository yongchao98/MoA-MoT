def find_r_flash_champions():
    """
    This function identifies champions from a predefined list who can perform
    the R-Flash combo in League of Legends up to Season 14.

    The R-Flash combo is defined as pressing R and then Flashing during the
    cast time or channel to reposition the ability's effect.
    """
    champions_who_can_r_flash = [
        "Amumu",
        "Braum",
        "Cassiopeia",
        "Diana",
        "Gnar",
        "Janna",
        "Lee Sin",
        "Neeko",
        "Nilah",
        "Poppy",
        "Qiyana",
        "Rell",
        "Riven",
        "Sejuani",
        "Seraphine",
        "Skarner",
        "Sona",
        "Vi",
        "Xin Zhao",
        "Yone"
    ]
    
    # The original list provided by the user for context, though not strictly needed for the output.
    all_champions_list = [
        "Amumu", "Aphelios", "Azir", "Blitzrank", "Braum", "Cassiopeia", 
        "Diana", "Evelynn", "Fizz", "Gnar", "Graves", "Gwen", "Hecarim", 
        "Illaoi", "Irelia", "Janna", "Jarvan", "Jax", "Kassadin", "Kindred", 
        "Lee Sin", "Malphite", "Milio", "Nami", "Nautilus", "Neeko", "Nilah", 
        "Oriana", "Poppy", "Qiyana", "Rell", "Renata", "Riven (for the 2nd R)", 
        "Sejuani", "Seraphine", "Skarner", "Sona", "Thresh", "Vi", "Xin Zhao", 
        "Yone", "Zehri"
    ]

    # For clarity, we'll correct "Riven (for the 2nd R)" to "Riven" as it is on our list.
    print("Champions who can R-Flash up until Season 14:")
    print(", ".join(champions_who_can_r_flash))

find_r_flash_champions()