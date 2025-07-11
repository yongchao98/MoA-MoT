def count_triple_crown_winners():
    """
    This script calculates and displays the total number of unique MLB players
    who have won the batting Triple Crown.
    """

    # Data for American League batting Triple Crown winners (Player, Year)
    al_winners = [
        ("Nap Lajoie", 1901),
        ("Ty Cobb", 1909),
        ("Jimmie Foxx", 1933),
        ("Lou Gehrig", 1934),
        ("Ted Williams", 1942),
        ("Ted Williams", 1947),  # Won twice
        ("Mickey Mantle", 1956),
        ("Frank Robinson", 1966),
        ("Carl Yastrzemski", 1967),
        ("Miguel Cabrera", 2012)
    ]

    # Data for National League batting Triple Crown winners (Player, Year)
    nl_winners = [
        ("Paul Hines", 1878),
        ("Tip O'Neill", 1887),
        ("Hugh Duffy", 1894),
        ("Heinie Zimmerman", 1912),
        ("Rogers Hornsby", 1922),
        ("Rogers Hornsby", 1925),  # Won twice
        ("Chuck Klein", 1933),
        ("Joe Medwick", 1937)
    ]

    # Get the unique player names from each league
    al_unique_players = set(player for player, year in al_winners)
    nl_unique_players = set(player for player, year in nl_winners)

    # Count the number of unique winners in each league
    num_al_winners = len(al_unique_players)
    num_nl_winners = len(nl_unique_players)

    # Total number of unique winners is the sum, as no player won in both leagues
    total_unique_winners = num_al_winners + num_nl_winners

    print("--- Unique MLB Batting Triple Crown Winners ---")
    
    print(f"\nAmerican League Winners ({num_al_winners} players):")
    for player in sorted(list(al_unique_players)):
        print(f"- {player}")

    print(f"\nNational League Winners ({num_nl_winners} players):")
    for player in sorted(list(nl_unique_players)):
        print(f"- {player}")

    print("\n--- Final Calculation ---")
    print("To find the total number of unique winners, we add the counts from each league.")
    print(f"Number of AL Winners + Number of NL Winners = Total Unique Winners")
    print(f"{num_al_winners} + {num_nl_winners} = {total_unique_winners}")


if __name__ == "__main__":
    count_triple_crown_winners()
