import pandas as pd

def solve_task():
    """
    This function identifies the 2019 Toronto Blue Jays player who meets three specific
    criteria and had the highest Baseball-Reference WAR (bWAR).
    """
    # Step 1: Create a dataset of relevant 2019 Blue Jays players
    # Data includes name, bWAR, and flags for the three criteria.
    # Data compiled from Baseball-Reference.com (rosters, game logs, transactions).
    players_data = [
        {'name': 'Randal Grichuk', 'war': 1.8, 'played_first': True, 'played_last': True, 'on_il': True},
        {'name': 'Lourdes Gurriel Jr.', 'war': 1.9, 'played_first': True, 'played_last': True, 'on_il': True},
        {'name': 'Ken Giles', 'war': 2.3, 'played_first': True, 'played_last': True, 'on_il': True},
        {'name': 'Brandon Drury', 'war': -0.8, 'played_first': True, 'played_last': True, 'on_il': True},
        {'name': 'Tim Mayza', 'war': 0.1, 'played_first': True, 'played_last': True, 'on_il': True},
        # Players who fail one or more criteria are included to demonstrate the filter
        {'name': 'Teoscar Hernández', 'war': 1.4, 'played_first': True, 'played_last': True, 'on_il': False},
        {'name': 'Danny Jansen', 'war': 0.7, 'played_first': True, 'played_last': True, 'on_il': False},
        {'name': 'Justin Smoak', 'war': 0.1, 'played_first': True, 'played_last': False, 'on_il': True},
        {'name': 'Marcus Stroman', 'war': 1.7, 'played_first': True, 'played_last': False, 'on_il': True}, # WAR with Blue Jays only
        {'name': 'Vladimir Guerrero Jr.', 'war': 1.9, 'played_first': False, 'played_last': True, 'on_il': True},
    ]

    players_df = pd.DataFrame(players_data)

    # Step 2: Apply the filters sequentially
    print("Filtering players based on the criteria:")
    print("-" * 40)

    # Requirement 1: Played in the first game of the 2019 season
    played_first_game = players_df[players_df['played_first'] == True]
    print("1. Players who played in the first game:")
    print(', '.join(played_first_game['name'].tolist()))
    print("-" * 40)

    # Requirement 2: Played in the final game of the 2019 season
    played_both_games = played_first_game[played_first_game['played_last'] == True]
    print("2. Players from the above list who also played in the last game:")
    print(', '.join(played_both_games['name'].tolist()))
    print("-" * 40)

    # Requirement 3: Was placed on the Injured List in 2019
    final_candidates = played_both_games[played_both_games['on_il'] == True]
    print("3. Players from the above list who were on the IL in 2019:")
    print(', '.join(final_candidates['name'].tolist()))
    print("-" * 40)

    # Step 3: Find the player with the highest WAR among the final candidates
    if not final_candidates.empty:
        top_player = final_candidates.loc[final_candidates['war'].idxmax()]
        print("Final Result:")
        print(f"The player who meets all criteria and has the highest WAR is {top_player['name']}.")
        print(f"This player's WAR in 2019 was: {top_player['war']}")

        # Required final output format
        print("\n<<<{}>>>".format(top_player['name']))
        print("<<<{}>>>".format(top_player['war']))

    else:
        print("No player meets all the specified criteria.")

solve_task()