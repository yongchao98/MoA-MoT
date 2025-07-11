import collections

def count_mlb_triple_crown_winners():
    """
    This script calculates and displays the number of unique MLB Triple Crown winners.
    It lists each unique player and shows a sum to represent the final count.
    """
    # Data for all MLB batting Triple Crown winning seasons.
    # Source: Major League Baseball official records.
    triple_crown_seasons = [
        {"player": "Paul Hines", "year": 1878},
        {"player": "Tip O'Neill", "year": 1887},
        {"player": "Hugh Duffy", "year": 1894},
        {"player": "Nap Lajoie", "year": 1901},
        {"player": "Ty Cobb", "year": 1909},
        {"player": "Rogers Hornsby", "year": 1922},
        {"player": "Rogers Hornsby", "year": 1925},
        {"player": "Chuck Klein", "year": 1933},
        {"player": "Jimmie Foxx", "year": 1933},
        {"player": "Lou Gehrig", "year": 1934},
        {"player": "Joe Medwick", "year": 1937},
        {"player": "Ted Williams", "year": 1942},
        {"player": "Ted Williams", "year": 1947},
        {"player": "Mickey Mantle", "year": 1956},
        {"player": "Frank Robinson", "year": 1966},
        {"player": "Carl Yastrzemski", "year": 1967},
        {"player": "Miguel Cabrera", "year": 2012},
    ]

    # Extract player names from the list of winning seasons
    player_names = [season['player'] for season in triple_crown_seasons]

    # A set is used to automatically handle duplicates, giving us a list of unique players
    unique_winners = sorted(list(set(player_names)))

    # The final answer is the count of unique players
    num_unique_winners = len(unique_winners)

    print("The unique players who have won the MLB Triple Crown are:")
    for name in unique_winners:
        print(f"- {name}")

    print("\nTo find the total number of unique winners, we sum them up.")
    print("The final equation representing the count of each unique player is:")
    
    # We represent each unique winner as '1' in our sum
    equation_parts = ['1' for _ in unique_winners]
    equation_str = " + ".join(equation_parts)

    print(f"{equation_str} = {num_unique_winners}")

    print(f"\nThere are {num_unique_winners} unique Triple Crown winners in MLB history.")

if __name__ == "__main__":
    count_mlb_triple_crown_winners()