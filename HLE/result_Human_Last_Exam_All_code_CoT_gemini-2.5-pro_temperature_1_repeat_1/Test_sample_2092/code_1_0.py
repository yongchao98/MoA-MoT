def find_r_flash_champions():
    """
    Identifies and prints the names of champions from a predefined list
    who can perform an R-Flash combo in League of Legends.
    """
    # The list of champions is based on an analysis of their ultimate abilities
    # and their interaction with the Flash summoner spell up to Season 14.
    # The R-Flash combo is defined as casting R then Flashing during the cast time
    # to make the ability originate from the new location.
    champion_list = [
        "Amumu",
        "Braum",
        "Cassiopeia",
        "Diana",
        "Evelynn",
        "Gnar",
        "Graves",
        "Hecarim",
        "Lee Sin",
        "Neeko",
        "Nilah",
        "Orianna",
        "Qiyana",
        "Rell",
        "Riven",
        "Sejuani",
        "Seraphine",
        "Skarner",
        "Sona",
        "Xin Zhao",
        "Yone"
    ]
    
    # Joining the list of names with ", " to create the final output string.
    result = ", ".join(champion_list)
    
    print(result)

find_r_flash_champions()
<<<Amumu, Braum, Cassiopeia, Diana, Evelynn, Gnar, Graves, Hecarim, Lee Sin, Neeko, Nilah, Orianna, Qiyana, Rell, Riven, Sejuani, Seraphine, Skarner, Sona, Xin Zhao, Yone>>>