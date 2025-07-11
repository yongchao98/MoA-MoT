def count_mlb_triple_crown_winners():
    """
    This function calculates and displays the total number of unique players
    who have won the MLB batting Triple Crown.
    """
    # Data representing each time the batting Triple Crown was won in the AL and NL.
    # Format: (Player Name, Year)
    al_wins = [
        ("Nap Lajoie", 1901), ("Ty Cobb", 1909), ("Jimmie Foxx", 1933),
        ("Lou Gehrig", 1934), ("Ted Williams", 1942), ("Ted Williams", 1947),
        ("Mickey Mantle", 1956), ("Frank Robinson", 1966),
        ("Carl Yastrzemski", 1967), ("Miguel Cabrera", 2012)
    ]

    nl_wins = [
        ("Paul Hines", 1878), ("Tip O'Neill", 1887), ("Hugh Duffy", 1894),
        ("Heinie Zimmerman", 1912), ("Rogers Hornsby", 1922),
        ("Rogers Hornsby", 1925), ("Chuck Klein", 1933), ("Joe Medwick", 1937)
    ]

    # Use a set to get the unique player names from each league.
    al_unique_winners = set(player for player, year in al_wins)
    nl_unique_winners = set(player for player, year in nl_wins)

    # Count the number of unique winners in each league.
    num_al_winners = len(al_unique_winners)
    num_nl_winners = len(nl_unique_winners)

    # Calculate the total number of unique winners.
    total_unique_winners = num_al_winners + num_nl_winners

    print("Calculating the number of unique MLB batting Triple Crown winners.")
    print("-" * 60)
    print(f"There are {num_al_winners} unique winners from the American League.")
    print(f"There are {num_nl_winners} unique winners from the National League.")
    print("-" * 60)
    print("The final equation to find the total number of winners is:")
    print(f"{num_al_winners} + {num_nl_winners} = {total_unique_winners}")
    print("-" * 60)
    print(f"In total, there have been {total_unique_winners} unique players who have won the MLB batting Triple Crown.")

count_mlb_triple_crown_winners()