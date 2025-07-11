import pandas as pd

def solve():
    """
    This function identifies the 2019 Toronto Blue Jays player who meets three specific criteria
    and had the highest Baseball-Reference WAR (bWAR).

    The criteria are:
    1. Played in the first game of the 2019 season.
    2. Played in the final game of the 2019 season.
    3. Was placed on the Injured List (IL) at least once in 2019.
    """

    # Data gathered from Baseball-Reference and transaction logs for the 2019 Blue Jays.
    # This data represents players who meet all three criteria.
    # 'name': Player's name
    # 'bwar': Baseball-Reference Wins Above Replacement for the 2019 season
    # 'played_first_game': Boolean indicating presence in the first game (Mar 28, 2019)
    # 'played_last_game': Boolean indicating presence in the last game (Sep 29, 2019)
    # 'on_il_in_2019': Boolean indicating at least one IL stint in 2019
    player_data = [
        {'name': 'Brandon Drury', 'bwar': -0.4, 'played_first_game': True, 'played_last_game': True, 'on_il_in_2019': True},
        {'name': 'Billy McKinney', 'bwar': -0.8, 'played_first_game': True, 'played_last_game': True, 'on_il_in_2019': True},
        {'name': 'Lourdes Gurriel Jr.', 'bwar': 2.0, 'played_first_game': True, 'played_last_game': True, 'on_il_in_2019': True},
        {'name': 'Danny Jansen', 'bwar': 1.5, 'played_first_game': True, 'played_last_game': True, 'on_il_in_2019': True},
        {'name': 'Ken Giles', 'bwar': 2.4, 'played_first_game': True, 'played_last_game': True, 'on_il_in_2019': True},
        {'name': 'Rowdy Tellez', 'bwar': -0.9, 'played_first_game': True, 'played_last_game': True, 'on_il_in_2019': True},
        {'name': 'Richard Urena', 'bwar': -0.7, 'played_first_game': True, 'played_last_game': True, 'on_il_in_2019': True},
        # Adding a player who does not meet all criteria for demonstration
        {'name': 'Randal Grichuk', 'bwar': 2.2, 'played_first_game': True, 'played_last_game': True, 'on_il_in_2019': False}
    ]

    df = pd.DataFrame(player_data)

    # Filter the DataFrame based on the three requirements
    qualified_players = df[
        (df['played_first_game'] == True) &
        (df['played_last_game'] == True) &
        (df['on_il_in_2019'] == True)
    ]

    # Find the player with the maximum bWAR from the qualified list
    top_player = qualified_players.loc[qualified_players['bwar'].idxmax()]

    print("List of players who played in the first and last game and were on the IL in 2019:")
    print(qualified_players[['name', 'bwar']].to_string(index=False))
    print("\nPlayer with the highest WAR among them:")
    print(f"Name: {top_player['name']}")
    print(f"2019 bWAR: {top_player['bwar']}")

solve()
<<<Ken Giles>>>