import pandas as pd

def solve():
    """
    This script identifies the 2019 Toronto Blue Jays player who meets three specific criteria
    and had the highest Baseball-Reference WAR (bWAR) for that season.
    """
    # Data manually gathered from Baseball-Reference.com for the 2019 season.

    # 1. Players in the first game of the season (March 28, 2019)
    first_game_players = {
        'Danny Jansen', 'Justin Smoak', 'Lourdes Gurriel Jr.', 'Freddy Galvis',
        'Brandon Drury', 'Randal Grichuk', 'Kevin Pillar', 'Teoscar Hernández',
        'Billy McKinney', 'Alen Hanson', 'Rowdy Tellez', 'Luke Maile',
        'Richard Ureña', 'Marcus Stroman', 'Joe Biagini', 'Tim Mayza',
        'Daniel Hudson', 'Ken Giles'
    }

    # 2. Players in the final game of the season (September 29, 2019)
    final_game_players = {
        'Reese McGuire', 'Rowdy Tellez', 'Cavan Biggio', 'Richard Ureña',
        'Brandon Drury', 'Billy McKinney', 'Jonathan Davis', 'Anthony Alford',
        'Breyvic Valera', 'Beau Taylor', 'Danny Jansen', 'Vladimir Guerrero Jr.',
        'Justin Smoak', 'Clay Buchholz', 'Sam Gaviglio', 'Buddy Boshers',
        'Brock Stewart', 'Jordan Romano', 'Ryan Tepera', 'Ken Giles'
    }

    # 3. Players placed on the Injured List during the 2019 season
    injured_list_players = {
        'Ryan Borucki', 'Matt Shoemaker', 'Clayton Richard', 'Devon Travis',
        'Dalton Pompey', 'David Phelps', 'John Axford', 'Bud Norris',
        'Julian Merryweather', 'Ryan Tepera', 'Elvis Luciano', 'Aaron Sanchez',
        'Marcus Stroman', 'Ken Giles', 'Justin Smoak', 'Danny Jansen',
        'Brandon Drury', 'Richard Ureña', 'Lourdes Gurriel Jr.', 'Luke Maile',
        'Bo Bichette', 'Jacob Waguespack', 'Tim Mayza'
    }
    
    # 4. Find the players who meet all three criteria using set intersection
    eligible_players = first_game_players.intersection(final_game_players).intersection(injured_list_players)
    
    print("Players who played in the first game, final game, AND were on the IL in 2019:")
    for player in sorted(list(eligible_players)):
        print(f"- {player}")
    print("-" * 30)

    # 5. WAR values (bWAR) for the eligible players
    player_war = {
        'Ken Giles': 2.4,
        'Danny Jansen': 1.0,
        'Justin Smoak': 0.3,
        'Brandon Drury': -0.9,
        'Richard Ureña': -0.5
    }

    # 6. Find the player with the highest WAR from the eligible list
    highest_war = -999
    top_player = None

    for player in eligible_players:
        if player in player_war:
            war = player_war[player]
            print(f"Checking {player}: WAR = {war}")
            if war > highest_war:
                highest_war = war
                top_player = player
    
    print("-" * 30)
    print(f"The player who meets all criteria with the highest WAR is {top_player}.")
    print(f"Their 2019 bWAR was: {highest_war}")


solve()
<<<Ken Giles>>>