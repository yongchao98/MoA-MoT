def count_triple_crown_winners():
    """
    This function calculates and displays the number of unique MLB batting Triple Crown winners.
    """
    # A list of tuples: (Player, Year) for each time the batting Triple Crown was won.
    triple_crown_instances = [
        ("Paul Hines", 1878),
        ("Tip O'Neill", 1887),
        ("Hugh Duffy", 1894),
        ("Nap Lajoie", 1901),
        ("Ty Cobb", 1909),
        ("Heinie Zimmerman", 1912),
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
        ("Miguel Cabrera", 2012),
    ]

    # Use a set to get the unique winners' names
    unique_winners = sorted(list(set(player for player, year in triple_crown_instances)))
    
    num_unique_winners = len(unique_winners)

    print("The unique MLB batting Triple Crown winners are:")
    
    # Create the equation string
    equation_parts = []
    for winner in unique_winners:
        # Each unique winner contributes 1 to the total count.
        equation_parts.append(f"1 ({winner})")
        
    equation_string = " + ".join(equation_parts)
    
    print(f"{equation_string} = {num_unique_winners}")
    print(f"\nIn total, there have been {num_unique_winners} unique Triple Crown winners in MLB history.")

count_triple_crown_winners()