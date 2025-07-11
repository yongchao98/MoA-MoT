# A list of unique players who have won the batting Triple Crown.
# A player is only listed once, even if they won multiple times (e.g., Ted Williams, Rogers Hornsby).
batting_winners = [
    "Paul Hines", "Tip O'Neill", "Hugh Duffy", "Nap Lajoie", "Ty Cobb",
    "Rogers Hornsby", "Chuck Klein", "Jimmie Foxx", "Lou Gehrig", "Joe Medwick",
    "Ted Williams", "Mickey Mantle", "Frank Robinson", "Carl Yastrzemski", "Miguel Cabrera"
]

# A list of unique players who have won the pitching Triple Crown.
# A player is only listed once, even if they won multiple times (e.g., Walter Johnson, Sandy Koufax).
pitching_winners = [
    "Tommy Bond", "Old Hoss Radbourn", "Guy Hecker", "Tim Keefe", "John Clarkson",
    "Amos Rusie", "Cy Young", "Rube Waddell", "Christy Mathewson", "Walter Johnson",
    "Grover Cleveland Alexander", "Dazzy Vance", "Lefty Grove", "Lefty Gomez", "Bucky Walters",
    "Bob Feller", "Hal Newhouser", "Sandy Koufax", "Steve Carlton", "Roger Clemens",
    "Randy Johnson", "Pedro Martinez", "Johan Santana", "Jake Peavy", "Clayton Kershaw",
    "Justin Verlander", "Shane Bieber", "Gerrit Cole"
]


num_batting_winners = len(batting_winners)
num_pitching_winners = len(pitching_winners)
total_winners = num_batting_winners + num_pitching_winners

print(f"Number of batting Triple Crown winners: {num_batting_winners}")
print(f"Number of pitching Triple Crown winners: {num_pitching_winners}")
print(f"Total MLB Triple Crown winners: {num_batting_winners} + {num_pitching_winners} = {total_winners}")