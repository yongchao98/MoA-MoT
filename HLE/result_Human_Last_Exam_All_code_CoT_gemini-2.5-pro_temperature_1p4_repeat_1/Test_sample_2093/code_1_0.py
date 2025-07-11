def solve_dash_flash():
    """
    This function identifies and prints the names of League of Legends champions
    from a predefined list who can perform the 'Dash-Flash' combo.
    The Dash-Flash combo is defined as using Flash during a dash ability's animation
    to buffer the ability and alter its landing spot and effect (damage/CC).
    """
    
    # List of champions who can perform the Dash-Flash combo as described.
    champions_with_combo = [
        "Gragas",
        "Vi",
        "Jarvan",
        "Galio",
        "Kayn",
        "Ornn",
        "Rakan",
        "Rell",
        "Riven",
        "Sejuani",
        "Shen",
        "Tryndamere",
        "Urgot",
        "Yone",
        "Zac",
        "Renekton",
        "Camille"
    ]
    
    # Join the list of names with a comma and space for clean output.
    result = ", ".join(champions_with_combo)
    
    # Print the final result.
    print(result)

solve_dash_flash()
