def find_dash_flash_champions():
    """
    This function identifies and prints the names of champions who can perform a "Dash-Flash" combo.
    The combo is defined as using Flash during a dash ability to buffer its effect (damage or CC)
    and apply it at the post-Flash location.
    """
    
    # List of champions who can perform the Dash-Flash combo as per the analysis
    champions = [
        "Gragas",
        "Vi",
        "Jarvan",
        "Galio",
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
        "Zac"
    ]
    
    # Join the list of names with a comma and print the result
    print(",".join(champions))

find_dash_flash_champions()