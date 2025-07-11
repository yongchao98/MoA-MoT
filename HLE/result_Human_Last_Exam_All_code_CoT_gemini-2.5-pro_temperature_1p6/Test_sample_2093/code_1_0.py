def solve():
    """
    This function identifies League of Legends champions from a predefined list 
    who can perform a 'Dash-Flash' combo. The criteria for this combo is that
    Flash is used during the dash ability to alter or extend the landing spot 
    of the ability's damage or crowd control, not just to reposition the champion.
    """
    
    # List of champions who can perform the Dash-Flash combo as defined.
    # Analysis summary:
    # Yes: Gragas, Vi, Jarvan, Galio, Kayn, Ornn, Rakan, Pyke, Rell, Riven, Sejuani, Shen, Tryndamere, Urgot, Yone, Renekton, Fiora, Camille, BelVeth.
    # No (repositioning only, targeted, or mechanically impossible): Fizz, LeBlanc, Poppy, Tristana, Zac, Nidalee, Corki.
    
    champions_who_can_dash_flash = [
        "Gragas",
        "Vi",
        "Jarvan",
        "Galio",
        "Kayn",
        "Ornn",
        "Rakan",
        "Pyke",
        "Rell",
        "Riven",
        "Sejuani",
        "Shen",
        "Tryndamere",
        "Urgot",
        "Yone",
        "Renekton",
        "Fiora",
        "Camille",
        "BelVeth"
    ]
    
    # Print the final list, with names separated by a comma and a space.
    result = ", ".join(champions_who_can_dash_flash)
    print(result)

solve()
<<<Gragas, Vi, Jarvan, Galio, Kayn, Ornn, Rakan, Pyke, Rell, Riven, Sejuani, Shen, Tryndamere, Urgot, Yone, Renekton, Fiora, Camille, BelVeth>>>