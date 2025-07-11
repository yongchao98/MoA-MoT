def count_triple_crown_winners():
    """
    This function calculates and prints the number of MLB Triple Crown winners.
    The Triple Crown is awarded to a batter who leads their league in
    batting average, home runs, and runs batted in (RBI) in the same season.
    """
    # Data of Triple Crown winners: Player -> List of years won
    # This structure correctly handles players who won multiple times.
    triple_crown_winners = {
        "Paul Hines": [1878],
        "Tip O'Neill": [1887],
        "Hugh Duffy": [1894],
        "Nap Lajoie": [1901],
        "Ty Cobb": [1909],
        "Heinie Zimmerman": [1912],
        "Rogers Hornsby": [1922, 1925],
        "Chuck Klein": [1933],
        "Jimmie Foxx": [1933],
        "Lou Gehrig": [1934],
        "Joe Medwick": [1937],
        "Ted Williams": [1942, 1947],
        "Mickey Mantle": [1956],
        "Frank Robinson": [1966],
        "Carl Yastrzemski": [1967],
        "Miguel Cabrera": [2012]
    }

    print("MLB Triple Crown Winners List:")
    print("-" * 30)
    
    # Sort the dictionary by the first year won for chronological order
    sorted_winners = sorted(triple_crown_winners.items(), key=lambda item: item[1][0])

    for player, years in sorted_winners:
        # Format years for display
        years_str = ", ".join(map(str, years))
        print(f"- {player} ({years_str})")

    # Calculate the total number of unique winners
    num_winners = len(triple_crown_winners)

    # Create the equation string "1 + 1 + ... + 1"
    equation = " + ".join(["1"] * num_winners)

    print("\nCalculating the total number of unique winners:")
    print(f"There are {equation} = {num_winners} unique Triple Crown winners in MLB history.")

count_triple_crown_winners()