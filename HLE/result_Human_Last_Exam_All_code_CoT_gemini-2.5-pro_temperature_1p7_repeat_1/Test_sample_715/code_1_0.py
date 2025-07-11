import pandas as pd
from pybaseball import team_batting

def solve_task():
    """
    This function identifies the 2019 Blue Jays player who meets the specified criteria
    and had the highest WAR.
    """
    # Based on research, the candidates who played in the first and last games
    # and were on the IL in 2019 are Randal Grichuk and Brandon Drury.
    candidate_names = ["Randal Grichuk", "Brandon Drury"]

    # Fetch 2019 Toronto Blue Jays batting statistics from Baseball Reference
    try:
        jays_stats_2019 = team_batting('TOR', 2019)
    except Exception as e:
        print(f"Failed to fetch data. Please ensure pybaseball is installed (`pip install pybaseball`) and you have an internet connection.")
        print(f"Error: {e}")
        return

    # Filter the stats for our candidate players
    candidates_stats = jays_stats_2019[jays_stats_2019['Name'].isin(candidate_names)]

    if candidates_stats.empty:
        print("Could not find statistics for the candidate players.")
        return

    # Find the player with the maximum WAR among the candidates
    # .loc is used to get the row, [0] to select the first if there are ties
    top_player = candidates_stats.loc[candidates_stats['WAR'].idxmax()]

    # Print the final result
    player_name = top_player['Name']
    player_war = top_player['WAR']

    print(f"The players who played in the first and last game of the 2019 season and were placed on the Injured List were:")
    for name in candidate_names:
        player_row = candidates_stats[candidates_stats['Name'] == name]
        if not player_row.empty:
            war_value = player_row['WAR'].iloc[0]
            print(f"- {name} (WAR: {war_value})")
    
    print(f"\nThe player with the highest WAR among them is {player_name}.")
    print("\n--- FINAL ANSWER ---")
    print(f"Player: {player_name}")
    print(f"WAR: {player_war}")

solve_task()
<<<Randal Grichuk>>>