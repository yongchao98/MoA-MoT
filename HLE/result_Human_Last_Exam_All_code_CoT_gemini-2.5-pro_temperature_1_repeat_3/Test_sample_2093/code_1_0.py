def solve():
    """
    Identifies and prints the names of League of Legends champions from a predefined list
    who can perform a 'Dash-Flash' combo.
    """
    # List of champions who can perform the Dash-Flash combo as described.
    # The list is ordered based on their appearance in the prompt.
    dash_flash_champions = [
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
        "Renekton",
        "Fiora"
    ]

    # Join the champion names into a single string, separated by ", "
    result = ", ".join(dash_flash_champions)
    
    # Print the final string
    print(result)

solve()
<<<Gragas, Vi, Jarvan, Galio, Ornn, Rakan, Pyke, Rell, Riven, Sejuani, Shen, Tryndamere, Urgot, Yone, Zac, Renekton, Fiora>>>