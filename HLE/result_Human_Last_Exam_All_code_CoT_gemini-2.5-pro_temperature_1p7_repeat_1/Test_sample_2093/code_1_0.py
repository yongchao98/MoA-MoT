def solve_dash_flash_champions():
    """
    Identifies and prints the names of League of Legends champions
    who can perform a buffered Dash-Flash combo as per the specified criteria.
    The criteria is that Flash is used during the dash animation to alter
    the ability's final impact zone, not just for repositioning after the fact.
    """

    # List of champions that can perform the defined Dash-Flash combo.
    # The ability is used, and during its travel/dash time, Flash can be
    # used to alter the final destination and effect of the ongoing ability.
    champions_with_dash_flash = [
        "Gragas",     # E-Flash: Can Flash during Body Slam to change stun location.
        "Vi",         # Q-Flash: Can Flash during Vault Breaker dash to redirect.
        "Jarvan",     # EQ-Flash: Can Flash during Dragon Strike dash to his flag.
        "Galio",      # E-Flash: Can Flash during Justice Punch dash.
        "Ornn",       # E-Flash: Can Flash during Searing Charge dash.
        "Rakan",      # W-Flash: Can Flash during Grand Entrance leap to change knock-up location.
        "Rell",       # W-Flash: Can Flash during Crash Down leap to change knock-up location.
        "Riven",      # Q3-Flash: Can Flash during the third cast of Broken Wings.
        "Sejuani",    # Q-Flash: Can Flash during Arctic Assault dash.
        "Shen",       # E-Flash: Can Flash during Shadow Dash to change taunt path/range.
        "Tryndamere", # E-Flash: Can Flash during Spinning Slash.
        "Urgot",      # E-Flash: Can Flash during Disdain charge to surprise-flip an enemy.
        "Zac",        # E-Flash: Can Flash mid-air during Elastic Slingshot.
        "Camille"     # E2-Flash: Can Flash during the Wall Dive portion of Hookshot.
    ]

    # Join the list into a single comma-separated string for the output.
    output_string = ",".join(champions_with_dash_flash)

    # Print the final result.
    print(output_string)

solve_dash_flash_champions()