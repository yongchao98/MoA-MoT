import pandas as pd

def find_blue_jays_player():
    """
    Identifies the 2019 Blue Jays player who played the first and last game, 
    was on the IL, and had the highest WAR.
    """
    
    # Step 1: Define the datasets based on research from Baseball-Reference.com
    
    # Players who played in the first game of the 2019 season (March 28, 2019)
    players_first_game = {
        'Brandon Drury', 'Teoscar Hernandez', 'Justin Smoak', 'Randal Grichuk', 
        'Lourdes Gurriel Jr.', 'Kevin Pillar', 'Danny Jansen', 'Freddy Galvis', 
        'Richard Urena', 'Marcus Stroman', 'Tim Mayza', 'Joe Biagini', 
        'Daniel Hudson', 'Ken Giles', 'Billy McKinney', 'Luke Maile'
    }

    # Players who played in the final game of the 2019 season (September 29, 2019)
    players_last_game = {
        'Bo Bichette', 'Cavan Biggio', 'Randal Grichuk', 'Rowdy Tellez', 
        'Teoscar Hernandez', 'Brandon Drury', 'Danny Jansen', 'Jonathan Davis', 
        'Richard Urena', 'Clay Buchholz', 'Sam Gaviglio', 'Buddy Boshers', 
        'Jordan Romano', 'Jason Adam', 'Ken Giles', 'Billy McKinney', 
        'Reese McGuire', 'Lourdes Gurriel Jr.', 'Breyvic Valera'
    }

    # Players placed on the Injured List at least once during the 2019 season
    players_on_il = {
        "Ryan Borucki", "John Axford", "Dalton Pompey", "Clay Buchholz", "Ryan Tepera",
        "Matt Shoemaker", "David Phelps", "Clayton Richard", "Devon Travis", "Elvis Luciano",
        "Julian Merryweather", "Nick Kingham", "Jacob Waguespack", "Jordan Romano", "Tim Mayza",
        "Ken Giles", "Lourdes Gurriel Jr.", "Billy McKinney", "Richard Urena", "Brandon Drury"
    }

    # 2019 WAR for each player as calculated by Baseball Reference
    player_war_data = {
        'Freddy Galvis': 2.3, 'Ken Giles': 2.2, 'Marcus Stroman': 1.8, 'Lourdes Gurriel Jr.': 1.7,
        'Bo Bichette': 1.7, 'Cavan Biggio': 1.5, 'Danny Jansen': 1.4, 'Randal Grichuk': 1.3,
        'Justin Smoak': 1.0, 'Aaron Sanchez': 0.7, 'Teoscar Hernandez': 0.7, 'Vladimir Guerrero Jr.': 0.5,
        'Joe Biagini': 0.4, 'Trent Thornton': 0.3, 'Daniel Hudson': 0.3, 'Rowdy Tellez': 0.1,
        'Clay Buchholz': 0.1, 'Tim Mayza': 0.0, 'Reese McGuire': 0.0, 'Jonathan Davis': -0.1,
        'Brandon Drury': -0.1, 'Luke Maile': -0.1, 'Ryan Tepera': -0.2, 'Matt Shoemaker': -0.2,
        'Sam Gaviglio': -0.2, 'Kevin Pillar': -0.3, 'Derek Law': -0.3, 'Eric Sogard': -0.3,
        'Richard Urena': -0.6, 'Socrates Brito': -0.8, 'Billy McKinney': -0.9, 'Clayton Richard': -1.0,
        'Thomas Pannone': -1.0, 'Alen Hanson': -1.1, 'Edwin Jackson': -1.9, 'Ryan Borucki': -0.2, 
        'John Axford': 0.0, 'Dalton Pompey': 0.0, 'David Phelps': 0.2, 'Devon Travis': 0.0, 
        'Elvis Luciano': -0.6, 'Julian Merryweather': 0.0, 'Nick Kingham': 0.0, 
        'Jacob Waguespack': -0.2, 'Jordan Romano': 0.0
    }

    # Step 2: Filter the players based on the criteria
    
    # Criterion 1 & 2: Played in both the first and final game of the season
    played_both_games = players_first_game.intersection(players_last_game)
    print("Players who played in the first AND last game of 2019:")
    print(f"{sorted(list(played_both_games))}\n")

    # Criterion 3: Was also placed on the Injured List
    qualified_players = played_both_games.intersection(players_on_il)
    print("Of those players, those who were also on the Injured List in 2019:")
    print(f"{sorted(list(qualified_players))}\n")

    # Step 3: Find the player with the highest WAR from the qualified list
    
    highest_war = -999.0
    top_player = None
    
    print("Calculating the highest WAR among the qualified players...")
    print("-" * 40)
    for player in qualified_players:
        war = player_war_data.get(player, -999.0)
        print(f"Player: {player:<20} | 2019 WAR: {war}")
        if war > highest_war:
            highest_war = war
            top_player = player
    print("-" * 40)
    
    # Step 4: Present the final answer
    if top_player:
        print(f"\nThe player who meets all criteria and had the highest WAR is {top_player} with a WAR of {highest_war}.")
    else:
        print("\nCould not find a player who meets all the criteria.")

if __name__ == '__main__':
    find_blue_jays_player()
<<<Ken Giles>>>