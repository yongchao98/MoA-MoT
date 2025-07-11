def find_dash_flash_champions():
    """
    Identifies and prints the names of champions from a predefined list who can
    perform the Dash-Flash combo.
    """
    # List of champions that can perform the Dash-Flash mechanic as described.
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
        "BelVeth"
    ]
    
    # Join the names with a comma and space for clean output.
    result = ", ".join(champions)
    
    print(result)

find_dash_flash_champions()