def find_dash_flash_champions():
    """
    This function identifies and prints the names of League of Legends champions
    from a predefined list who can perform the "Dash-Flash" combo.
    """
    # List of champions who can perform the Dash-Flash combo as described.
    champions = [
        "Gragas",
        "Vi",
        "Jarvan",
        "Galio",
        "Ornn",
        "Rakan",
        "Rell",
        "Riven",
        "Sejuani",
        "Shen",
        "Tryndamere",
        "Urgot",
        "Yone",
        "Camille"
    ]

    # Join the list elements into a single string, separated by commas.
    result = ",".join(champions)

    # Print the final result.
    print(result)

find_dash_flash_champions()