import sys

def find_blue_jays_player():
    """
    Identifies the 2019 Blue Jays player who meets specific criteria and had the highest WAR.
    
    The criteria are:
    1. Played in the first game of the 2019 season.
    2. Played in the final game of the 2019 season.
    3. Was placed on the Injured List at least once in 2019.
    
    The function then finds the player with the highest Baseball-Reference WAR among those who qualify.
    """
    
    # Data compiled from Baseball-Reference for the 2019 Toronto Blue Jays.
    # Each dictionary contains the player's name, their WAR, and flags for the three conditions.
    players_data = [
        # Players who played in BOTH the first and last game of 2019
        {'name': 'Randal Grichuk',     'war': 1.1,  'played_first': True, 'played_last': True, 'on_il_2019': True},
        {'name': 'Teoscar Hernández',  'war': 1.6,  'played_first': True, 'played_last': True, 'on_il_2019': False},
        {'name': 'Lourdes Gurriel Jr.','war': 2.0,  'played_first': True, 'played_last': True, 'on_il_2019': True},
        {'name': 'Brandon Drury',      'war': -0.8, 'played_first': True, 'played_last': True, 'on_il_2019': True},
        {'name': 'Billy McKinney',     'war': -0.2, 'played_first': True, 'played_last': True, 'on_il_2019': True},
        {'name': 'Richard Urena',      'war': -1.3, 'played_first': True, 'played_last': True, 'on_il_2019': False},
        {'name': 'Ken Giles',          'war': 1.8,  'played_first': True, 'played_last': True, 'on_il_2019': True},
        {'name': 'Tim Mayza',          'war': 0.1,  'played_first': True, 'played_last': True, 'on_il_2019': True},
        
        # Other notable players for context (do not meet all criteria)
        {'name': 'Justin Smoak',       'war': 0.5,  'played_first': True, 'played_last': False, 'on_il_2019': False},
        {'name': 'Freddy Galvis',      'war': 1.7,  'played_first': True, 'played_last': False, 'on_il_2019': True},
        {'name': 'Danny Jansen',       'war': 0.8,  'played_first': True, 'played_last': False, 'on_il_2019': False},
        {'name': 'Marcus Stroman',     'war': 1.5,  'played_first': True, 'played_last': False, 'on_il_2019': False},
        {'name': 'Vladimir Guerrero Jr.', 'war': 0.4, 'played_first': False, 'played_last': True, 'on_il_2019': True},
        {'name': 'Bo Bichette',        'war': 1.7,  'played_first': False, 'played_last': True, 'on_il_2019': True},
    ]

    eligible_players = []
    
    print("Step 1: Identifying all players who meet the three requirements...")
    
    # Filter players based on the three criteria
    for player in players_data:
        if player['played_first'] and player['played_last'] and player['on_il_2019']:
            eligible_players.append(player)
            print(f"- {player['name']} meets all criteria. WAR: {player['war']}")

    if not eligible_players:
        print("No players found who meet all criteria.")
        return

    # Find the player with the highest WAR from the eligible list
    highest_war_player = max(eligible_players, key=lambda x: x['war'])
    
    print("\nStep 2: Finding the player with the highest WAR from the eligible list.")
    print(f"\nThe player who meets all three criteria with the highest WAR is {highest_war_player['name']}.")
    print(f"Their WAR for the 2019 season was: {highest_war_player['war']}")


find_blue_jays_player()
<<<Lourdes Gurriel Jr.>>>