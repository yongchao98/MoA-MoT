def find_blue_jays_player():
    """
    Identifies the 2019 Blue Jays player who meets specific criteria
    and had the highest WAR.
    """
    # Step 1: Define the player data for the 2019 Toronto Blue Jays.
    # Data includes name, if they played in the first game, last game,
    # were on the Injured List (IL), and their 2019 WAR from Baseball Reference.
    players_data = [
        {'name': 'Brandon Drury', 'first_game': True, 'last_game': True, 'on_il': True, 'war': -0.4},
        {'name': 'Billy McKinney', 'first_game': True, 'last_game': True, 'on_il': True, 'war': -0.8},
        {'name': 'Teoscar Hernández', 'first_game': True, 'last_game': True, 'on_il': False, 'war': 1.6},
        {'name': 'Justin Smoak', 'first_game': True, 'last_game': True, 'on_il': True, 'war': 0.2},
        {'name': 'Randal Grichuk', 'first_game': True, 'last_game': True, 'on_il': False, 'war': 1.5},
        {'name': 'Danny Jansen', 'first_game': True, 'last_game': True, 'on_il': False, 'war': 1.4},
        {'name': 'Ken Giles', 'first_game': True, 'last_game': True, 'on_il': True, 'war': 2.3},
        {'name': 'Rowdy Tellez', 'first_game': True, 'last_game': True, 'on_il': True, 'war': -1.1},
        {'name': 'Richard Ureña', 'first_game': True, 'last_game': True, 'on_il': True, 'war': -0.6},
        # Other players listed for context who do not meet all criteria
        {'name': 'Lourdes Gurriel Jr.', 'first_game': True, 'last_game': False, 'on_il': True, 'war': 1.7},
        {'name': 'Freddy Galvis', 'first_game': True, 'last_game': False, 'on_il': True, 'war': 1.8},
        {'name': 'Marcus Stroman', 'first_game': True, 'last_game': False, 'on_il': False, 'war': 3.1},
        {'name': 'Bo Bichette', 'first_game': False, 'last_game': True, 'on_il': True, 'war': 1.7},
        {'name': 'Cavan Biggio', 'first_game': False, 'last_game': True, 'on_il': True, 'war': 2.0},
    ]

    # Step 2: Filter for players who meet all three requirements.
    eligible_players = [
        player for player in players_data
        if player['first_game'] and player['last_game'] and player['on_il']
    ]

    # Step 3: Identify the player with the highest WAR from the eligible list.
    if not eligible_players:
        print("No player was found who meets all the specified criteria.")
        return

    print("The following players meet all three requirements:")
    print("1) Played in the first game of the 2019 season.")
    print("2) Played in the final game of the 2019 season.")
    print("3) Was placed on the Injured List in 2019.")
    print("-" * 30)

    # Sort players by WAR in descending order to create a ranking
    eligible_players.sort(key=lambda p: p['war'], reverse=True)

    # Print each eligible player and their WAR
    for player in eligible_players:
        print(f"Player: {player['name']}, WAR: {player['war']}")

    # The player with the highest WAR is the first one in the sorted list.
    top_player = eligible_players[0]

    print("\n" + "=" * 30)
    print(f"The player with the most WAR among this group is {top_player['name']}.")
    print(f"Final Answer: {top_player['name']} accumulated a WAR of {top_player['war']}.")
    print("=" * 30)

if __name__ == '__main__':
    find_blue_jays_player()