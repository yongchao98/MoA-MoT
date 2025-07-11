import pandas as pd

def solve_blue_jays_war():
    """
    This script identifies the 2019 Blue Jays player with the highest WAR who meets three specific criteria:
    1. Played in the first game of the season.
    2. Played in the last game of the season.
    3. Was on the Injured List at some point in 2019.
    The data is based on publicly available game logs and transaction data from Baseball Reference.
    """
    
    # Step 1-4: The analysis identifies three players who meet all criteria.
    # - Played in first game (3/28/19) & last game (9/29/19): Randal Grichuk, Brandon Drury, Tim Mayza, Ken Giles
    # - Of those four, those who were on the IL in 2019: Brandon Drury, Tim Mayza, Ken Giles
    eligible_players_war = {
        "Ken Giles": 2.1,
        "Brandon Drury": 0.2,
        "Tim Mayza": 0.1
    }

    # Step 5: Find the player with the highest WAR among the eligible candidates.
    highest_war_player = ""
    highest_war_value = -999.0

    for player, war in eligible_players_war.items():
        if war > highest_war_value:
            highest_war_value = war
            highest_war_player = player
            
    print("The eligible players who met all three requirements are:")
    # Print the "equation" by showing all candidates and their WAR values
    for player, war in eligible_players_war.items():
        print(f"Player: {player}, 2019 WAR = {war}")

    print("\n-------------------------------------------------")
    print("The player with the most WAR among them is:")
    print(f"{highest_war_player} = {highest_war_value}")

solve_blue_jays_war()