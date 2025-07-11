def solve_r_flash_champions():
    """
    This function identifies and prints the names of League of Legends champions
    from a predefined list who can perform an R-Flash combo up until Season 14.
    An R-Flash is defined as buffering the Flash spell during the ultimate ability's
    cast time or channel to reposition its effect.
    """

    # List of all champions provided in the prompt.
    all_champions = [
        "Amumu", "Aphelios", "Azir", "Blitzrank", "Braum", "Cassiopeia", "Diana",
        "Evelynn", "Fizz", "Gnar", "Graves", "Gwen", "Hecarim", "Illaoi", "Irelia",
        "Janna", "Jarvan", "Jax", "Kassadin", "Kindred", "Lee Sin", "Malphite",
        "Milio", "Nami", "Nautilus", "Neeko", "Nilah", "Oriana", "Poppy",
        "Qiyana", "Rell", "Renata", "Riven (for the 2nd R)", "Sejuani", "Seraphine",
        "Skarner", "Sona", "Thresh", "Vi", "Xin Zhao", "Yone", "Zehri"
    ]

    # Champions who can perform the R-Flash combo as defined.
    # The reasoning for each is based on their R's mechanics:
    # - Amumu, Braum, Cassiopeia, Gnar, Lee Sin, Neeko, Qiyana, Riven, Sejuani, Seraphine, Sona, Thresh, Yone:
    #   Can Flash during the R cast/channel animation to change its origin point or angle.
    # - Azir: Can R-Flash to change the wall's starting position (Shurima Shuffle).
    # - Diana, Milio, Nilah, Rell, Xin Zhao: Can R-Flash to move the center of their radial AOE pull/push/heal.
    # - Skarner (pre-rework): Can R a target and Flash to drag them further.
    r_flash_champions = [
        "Amumu",
        "Azir",
        "Braum",
        "Cassiopeia",
        "Diana",
        "Gnar",
        "Lee Sin",
        "Milio",
        "Neeko",
        "Nilah",
        "Qiyana",
        "Rell",
        "Riven",
        "Sejuani",
        "Seraphine",
        "Skarner",
        "Sona",
        "Thresh",
        "Xin Zhao",
        "Yone",
    ]

    # In the problem, "Riven (for the 2nd R)" is a special case. We will use "Riven" as the name.
    
    # We will also use this list to generate the final output string.
    final_champions_list = sorted(r_flash_champions)

    print(",".join(final_champions_list))

solve_r_flash_champions()
<<<Amumu,Azir,Braum,Cassiopeia,Diana,Gnar,Lee Sin,Milio,Neeko,Nilah,Qiyana,Rell,Riven,Sejuani,Seraphine,Skarner,Sona,Thresh,Xin Zhao,Yone>>>