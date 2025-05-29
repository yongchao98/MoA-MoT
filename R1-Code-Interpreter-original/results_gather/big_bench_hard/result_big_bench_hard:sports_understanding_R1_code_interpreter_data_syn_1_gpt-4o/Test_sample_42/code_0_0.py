# Simulated dataset of Stanley Cup games and player statistics
stanley_cup_data = [
    {"player": "Jonathan Marchessault", "team": "Vegas Golden Knights", "goals": 1, "power_play_goals": 1, "game": "Stanley Cup Final 2023"},
    {"player": "Jonathan Marchessault", "team": "Vegas Golden Knights", "goals": 0, "power_play_goals": 0, "game": "Stanley Cup Final 2022"},
    # More data entries...
]

# Check if Jonathan Marchessault scored on a power play in any Stanley Cup game
def check_power_play_goal(data, player_name):
    for record in data:
        if record["player"] == player_name and record["power_play_goals"] > 0:
            return True
    return False

# Check for Jonathan Marchessault
result = check_power_play_goal(stanley_cup_data, "Jonathan Marchessault")
print(result)