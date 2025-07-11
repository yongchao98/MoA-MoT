import pandas as pd
import io

def solve():
    """
    This function identifies the 2019 Toronto Blue Jays player who meets the specified criteria
    and had the highest WAR as calculated by Baseball Reference.
    """
    # Criterion 1: Players in the first game of the 2019 season (March 28, 2019)
    first_game_players = [
        "Billy McKinney", "Teoscar Hernández", "Randal Grichuk", "Justin Smoak",
        "Kevin Pillar", "Lourdes Gurriel Jr.", "Brandon Drury", "Danny Jansen",
        "Freddy Galvis", "Marcus Stroman", "Tim Mayza", "Joe Biagini", "Ken Giles", "Luke Maile"
    ]

    # Criterion 2: Players in the final game of the 2019 season (September 29, 2019)
    final_game_players = [
        "Bo Bichette", "Cavan Biggio", "Randal Grichuk", "Rowdy Tellez",
        "Teoscar Hernández", "Brandon Drury", "Jonathan Davis", "Danny Jansen",
        "Richard Ureña", "Clay Buchholz", "T.J. Zeuch", "Buddy Boshers",
        "Jason Adam", "Ken Giles", "Jordan Romano", "Justin Shafer", "Breyvic Valera"
    ]

    # Criterion 3: Players placed on the Injured List during the 2019 season
    # This list is based on transaction data for the 2019 Blue Jays.
    injured_list_players = [
        "Ryan Tepera", "Clayton Richard", "Devon Travis", "John Axford",
        "Bud Norris", "Dalton Pompey", "Ryan Borucki", "Matt Shoemaker",
        "David Phelps", "Aaron Sanchez", "Brandon Drury", "Randal Grichuk",
        "Elvis Luciano", "Lourdes Gurriel Jr.", "Tim Mayza", "Ken Giles",
        "Luke Maile", "Bo Bichette", "Jordan Romano"
    ]

    # Find the intersection of all three lists
    eligible_players = list(
        set(first_game_players) & set(final_game_players) & set(injured_list_players)
    )

    # Player WAR data from Baseball-Reference for the 2019 season
    # Using a dictionary for easy lookup: {Player Name: WAR}
    player_war = {
        "Randal Grichuk": 1.0,
        "Brandon Drury": -0.7,
        "Ken Giles": 2.3
        # Other players' WAR is not needed as they don't meet all criteria
    }

    # Find the eligible player with the highest WAR
    highest_war = -999
    top_player = None

    for player in eligible_players:
        if player in player_war and player_war[player] > highest_war:
            highest_war = player_war[player]
            top_player = player

    if top_player:
        print(f"The player who meets all three requirements is {top_player}.")
        print(f"His WAR for the 2019 season was {highest_war}.")
    else:
        print("No player was found who met all three requirements.")

solve()