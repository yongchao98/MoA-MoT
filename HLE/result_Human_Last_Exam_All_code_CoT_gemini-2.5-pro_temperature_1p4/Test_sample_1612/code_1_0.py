def find_triple_crown_winners():
    """
    This script lists the MLB Triple Crown winners and calculates the total number.
    The data includes winners from the American League, National League, and Negro Leagues,
    whose statistics are now officially recognized by MLB.
    """
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
        {"player": "Aaron Judge", "year": 2022}
    ]

    nl_winners = [
        {"player": "Heinie Zimmerman", "year": 1912},
        {"player": "Rogers Hornsby", "year": 1922},
        {"player": "Rogers Hornsby", "year": 1925},
        {"player": "Chuck Klein", "year": 1933},
        {"player": "Joe Medwick", "year": 1937}
    ]
    
    nel_winners = [
        {"player": "Oscar Charleston", "year": 1921},
        {"player": "Turkey Stearnes", "year": 1928},
        {"player": "Josh Gibson", "year": 1936},
        {"player": "Josh Gibson", "year": 1937},
        {"player": "Ted Strong", "year": 1939},
        {"player": "Lennie Pearson", "year": 1940},
        {"player": "Wild Bill Wright", "year": 1940},
        {"player": "Artie Wilson", "year": 1948}
    ]

    print("--- MLB Triple Crown Winners ---")
    
    print("\nAmerican League Winners:")
    for winner in al_winners:
        print(f"- {winner['player']} ({winner['year']})")

    print("\nNational League Winners:")
    for winner in nl_winners:
        print(f"- {winner['player']} ({winner['year']})")
        
    print("\nNegro Leagues Winners:")
    for winner in nel_winners:
        print(f"- {winner['player']} ({winner['year']})")

    num_al_winners = len(al_winners)
    num_nl_winners = len(nl_winners)
    num_nel_winners = len(nel_winners)
    total_winners = num_al_winners + num_nl_winners + num_nel_winners

    print("\n--- Total Number of Triple Crown Winners ---")
    print(f"American League Wins: {num_al_winners}")
    print(f"National League Wins: {num_nl_winners}")
    print(f"Negro Leagues Wins: {num_nel_winners}")
    print(f"\nFinal Equation: {num_al_winners} + {num_nl_winners} + {num_nel_winners} = {total_winners}")

find_triple_crown_winners()
<<<24>>>