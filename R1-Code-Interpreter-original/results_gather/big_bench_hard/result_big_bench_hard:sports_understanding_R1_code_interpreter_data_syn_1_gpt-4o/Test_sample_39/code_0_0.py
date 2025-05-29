import requests

def check_player(name):
    # Using a simple API or database lookup to check if Robert Woods is a known football player
    # For demonstration, we'll use a mock check
    known_players = ["Robert Woods"]
    return name in known_players

player_name = "Robert Woods"
is_player = check_player(player_name)
print(is_player)