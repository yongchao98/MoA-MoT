def find_dash_flash_champions():
    """
    This function identifies and prints the names of League of Legends champions
    from a predefined list who can perform a specific 'Dash-Flash' combo.
    
    The combo is defined as using Flash during a dash ability's animation
    to alter its landing point and more effectively apply the ability's own
    damage or crowd control.
    """
    
    # Based on in-game mechanics analysis, the following champions can perform
    # the described buffered Dash-Flash combo with the specified ability.
    champions = [
        "Gragas",      # E - Body Slam
        "Vi",          # Q - Vault Breaker
        "Jarvan",      # E+Q - Dragon Strike
        "Galio",       # E - Justice Punch
        "Ornn",        # E - Searing Charge
        "Rakan",       # W - Grand Entrance
        "Rell",        # W - Ferromancy: Crash Down
        "Riven",       # Q3 - Broken Wings
        "Sejuani",     # Q - Arctic Assault
        "Shen",        # E - Shadow Dash
        "Tryndamere",  # E - Spinning Slash
        "Urgot",       # E - Disdain
        "Fiora",       # Q - Lunge
        "Camille",     # E - Wall Dive
        "BelVeth"      # Q - Void Surge
    ]
    
    # Join the list of champion names with a comma and print the result.
    print(",".join(champions))

find_dash_flash_champions()