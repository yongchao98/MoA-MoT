def count_mlb_triple_crown_winners():
    """
    This function identifies and counts the unique players who have won
    the MLB Triple Crown for batting.
    """
    # A list of tuples, where each tuple is (Player, Year of Win)
    triple_crown_wins = [
        ("Paul Hines", 1878),
        ("Tip O'Neill", 1887),
        ("Hugh Duffy", 1894),
        ("Nap Lajoie", 1901),
        ("Ty Cobb", 1909),
        ("Rogers Hornsby", 1922),
        ("Rogers Hornsby", 1925),
        ("Chuck Klein", 1933),
        ("Jimmie Foxx", 1933),
        ("Lou Gehrig", 1934),
        ("Joe Medwick", 1937),
        ("Ted Williams", 1942),
        ("Ted Williams", 1947),
        ("Mickey Mantle", 1956),
        ("Frank Robinson", 1966),
        ("Carl Yastrzemski", 1967),
        ("Miguel Cabrera", 2012)
    ]

    # Use a set to get a list of unique player names, then sort it
    unique_winners = sorted(list(set(player for player, year in triple_crown_wins)))
    
    print("The unique MLB Triple Crown winners are:")
    for name in unique_winners:
        print(name)
        
    count = len(unique_winners)
    
    # Create the equation string as requested
    equation_str = " + ".join(['1'] * count)

    print("\nThe total number of unique winners is found by summing the count of each player:")
    # We output each number in the final equation, which are all '1's representing each unique player.
    print(f"{equation_str} = {count}")

count_mlb_triple_crown_winners()