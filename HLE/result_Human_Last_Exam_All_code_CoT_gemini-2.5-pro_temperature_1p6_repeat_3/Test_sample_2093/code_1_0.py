def find_dash_flash_champs():
    """
    Identifies champions from a predefined list who can perform a Dash-Flash combo
    where Flash is used to buffer and redirect the dash's primary effect (damage/CC).
    """
    
    # List of champions that can perform the specified Dash-Flash combo
    capable_champions = [
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
        "Zac",
        "Camille"
    ]
    
    # Join the list into a single string with commas
    result = ",".join(capable_champions)
    
    # Print the final result
    print(result)

find_dash_flash_champs()