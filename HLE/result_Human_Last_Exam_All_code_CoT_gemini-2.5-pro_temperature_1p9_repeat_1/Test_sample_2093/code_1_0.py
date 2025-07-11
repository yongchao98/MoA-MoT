def find_dash_flash_champions():
    """
    This function identifies and lists League of Legends champions
    who can perform a buffered Dash-Flash combo with the specified ability.
    """
    
    # List of champions who can perform the buffered Dash-Flash combo
    # based on the provided list and ability.
    champions_with_combo = [
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
        "Zac",
        "Renekton",
        "Fiora",
        "Camille",
        "BelVeth"
    ]
    
    # Sort the list alphabetically for clean output
    champions_with_combo.sort()

    # Join the names with a comma as the separator
    result = ",".join(champions_with_combo)
    
    # Print the final result
    print(result)

find_dash_flash_champions()