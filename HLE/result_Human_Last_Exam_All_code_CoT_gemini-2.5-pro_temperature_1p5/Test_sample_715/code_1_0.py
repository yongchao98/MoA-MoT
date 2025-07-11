import pandas as pd

def solve_blue_jays_riddle():
    """
    This script identifies a 2019 Toronto Blue Jays player based on three criteria:
    1. Played in the first game of the 2019 season.
    2. Played in the final game of the 2019 season.
    3. Was placed on the Injured List at least once during the 2019 season.
    
    It then finds the player who met these criteria and had the highest WAR
    as calculated by Baseball-Reference.com.
    """

    # Data gathered from Baseball-Reference.com
    # Columns: Player Name, Played in First Game, Played in Last Game, On Injured List, 2019 bWAR
    data = {
        'player': [
            'Randal Grichuk', 'Rowdy Tellez', 'Brandon Drury', 'Danny Jansen',
            'Ken Giles', 'Billy McKinney', 'Richard Urena', 'Justin Smoak',
            'Freddy Galvis', 'Bo Bichette', 'Cavan Biggio'
        ],
        'played_first_game': [
            True, True, True, True, True, True, True, True, True, False, False
        ],
        'played_last_game': [
            True, True, True, True, True, True, True, False, False, True, True
        ],
        'was_on_il': [
            True, False, True, True, True, True, True, True, True, True, True
        ],
        'bwar_2019': [
            1.5, 0.0, -0.3, 1.0, 1.9, -0.8, -0.3, 0.1, 1.6, 1.7, 2.2
        ]
    }
    
    df = pd.DataFrame(data)

    # Apply the three requirements to filter the DataFrame
    qualifying_players = df[
        (df['played_first_game'] == True) &
        (df['played_last_game'] == True) &
        (df['was_on_il'] == True)
    ]

    if qualifying_players.empty:
        print("No players were found who meet all three criteria.")
        return

    # Find the player with the maximum WAR among the qualified players
    # The idxmax() function returns the index of the first occurrence of the maximum value
    best_player_index = qualifying_players['bwar_2019'].idxmax()
    best_player = df.loc[best_player_index]

    player_name = best_player['player']
    player_war = best_player['bwar_2019']

    print(f"The qualifying player with the highest WAR in 2019 is {player_name}.")
    print(f"Their Baseball-Reference WAR was: {player_war}")

solve_blue_jays_riddle()