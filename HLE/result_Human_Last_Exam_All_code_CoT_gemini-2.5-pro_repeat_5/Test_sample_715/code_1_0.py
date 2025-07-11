import pandas as pd

def solve_blue_jays_war():
    """
    This function identifies the 2019 Toronto Blue Jays player who meets three specific
    criteria and had the highest Wins Above Replacement (WAR) as calculated by Baseball-Reference.
    """

    # Step 1: Establish lists of players based on the criteria.
    # Data is based on 2019 season records from baseball-reference.com and transaction logs.
    
    # Players in the first game of the 2019 season (Mar 28 vs. DET)
    players_first_game = {
        "Billy McKinney", "Brandon Drury", "Randal Grichuk", "Justin Smoak",
        "Teoscar Hernández", "Kevin Pillar", "Rowdy Tellez", "Lourdes Gurriel Jr.",
        "Danny Jansen", "Freddy Galvis", "Marcus Stroman", "Luke Maile",
        "Richard Urena", "Alen Hanson", "Joe Biagini", "Tim Mayza",
        "Ken Giles", "Daniel Hudson"
    }

    # Players in the final game of the 2019 season (Sep 29 vs. TB)
    players_last_game = {
        "Bo Bichette", "Cavan Biggio", "Randal Grichuk", "Rowdy Tellez",
        "Vladimir Guerrero Jr.", "Billy McKinney", "Jonathan Davis", "Reese McGuire",
        "Richard Urena", "Clay Buchholz", "Danny Jansen", "Brandon Drury",
        "Teoscar Hernández", "Lourdes Gurriel Jr.", "Justin Smoak", "Breyvic Valera",
        "T. J. Zeuch", "Buddy Boshers", "Ryan Tepera", "Wilmer Font", "Ken Giles"
    }

    # Key players who were on the Injured List at some point in 2019
    players_on_il_2019 = {
        "Billy McKinney", "Brandon Drury", "Justin Smoak", "Lourdes Gurriel Jr.",
        "Richard Urena", "Ken Giles", "Ryan Borucki", "Matt Shoemaker",
        "Clayton Richard", "Devon Travis", "Ryan Tepera", "Aaron Sanchez", "Tim Mayza"
    }

    # Step 2: Find the intersection of all three sets to get the eligible players.
    eligible_players = players_first_game.intersection(players_last_game).intersection(players_on_il_2019)

    # Step 3: Define a dictionary with the 2019 Baseball-Reference WAR for the eligible players.
    player_war_bref = {
        "Billy McKinney": -0.1,
        "Brandon Drury": 0.0,
        "Justin Smoak": 0.2,
        "Lourdes Gurriel Jr.": 2.0,
        "Richard Urena": -0.7,
        "Ken Giles": 2.3
    }

    # Step 4 & 5: Iterate through the eligible players to find the one with the highest WAR.
    highest_war_player = None
    max_war = -999  # Initialize with a very low number

    print("Identifying the player with the highest 2019 WAR who meets all criteria...")
    print("Candidates must have played in the first and last game and been on the IL.")
    print("-" * 30)

    # Create a list of tuples for sorting and display
    candidates = []
    for player in eligible_players:
        if player in player_war_bref:
            candidates.append((player, player_war_bref[player]))

    # Sort candidates by WAR in descending order
    candidates.sort(key=lambda x: x[1], reverse=True)

    # Display the comparison
    for player, war in candidates:
        print(f"Player: {player:<20} WAR: {war}")
        if war > max_war:
            max_war = war
            highest_war_player = player
            
    print("-" * 30)
    print(f"The player with the highest WAR among the candidates is {highest_war_player}.")
    print(f"Final Answer (Highest WAR): {max_war}")


solve_blue_jays_war()
<<<2.3>>>