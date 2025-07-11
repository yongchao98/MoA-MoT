def count_triple_crown_winners():
    """
    Calculates and displays the number of MLB Triple Crown winners.
    The Triple Crown is awarded to a batter who leads their league in
    batting average, home runs, and runs batted in (RBI) in the same season.
    """

    # Data for Triple Crown winners in the American League
    al_winners = [
        {"player": "Nap Lajoie", "year": 1901},
        {"player": "Ty Cobb", "year": 1909},
        {"player": "Jimmie Foxx", "year": 1933},
        {"player": "Lou Gehrig", "year": 1934},
        {"player": "Ted Williams", "year": 1942},
        {"player": "Ted Williams", "year": 1947},
        {"player": "Mickey Mantle", "year": 1956},
        {"player": "Frank Robinson", "year": 1966},
        {"player": "Carl Yastrzemski", "year": 1967},
        {"player": "Miguel Cabrera", "year": 2012},
    ]

    # Data for Triple Crown winners in the National League (and predecessor leagues)
    nl_winners = [
        {"player": "Paul Hines", "year": 1878},
        {"player": "Tip O'Neill", "year": 1887},
        {"player": "Hugh Duffy", "year": 1894},
        {"player": "Rogers Hornsby", "year": 1922},
        {"player": "Rogers Hornsby", "year": 1925},
        {"player": "Chuck Klein", "year": 1933},
        {"player": "Joe Medwick", "year": 1937},
    ]

    num_al_winners = len(al_winners)
    num_nl_winners = len(nl_winners)
    total_winners = num_al_winners + num_nl_winners

    print(f"There have been {num_al_winners} Triple Crown winning seasons in the American League.")
    print(f"There have been {num_nl_winners} Triple Crown winning seasons in the National League.")
    print("\nTo find the total number of Triple Crown winners, we add the winners from both leagues:")
    print(f"Final Equation: {num_al_winners} + {num_nl_winners} = {total_winners}")

count_triple_crown_winners()