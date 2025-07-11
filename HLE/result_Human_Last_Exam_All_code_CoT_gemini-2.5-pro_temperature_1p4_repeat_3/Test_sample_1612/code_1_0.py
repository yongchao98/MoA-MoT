def count_triple_crown_winners():
    """
    Calculates and prints the number of MLB Triple Crown winners.
    The function stores the data of all Triple Crown instances,
    calculates the totals per league, the total instances, and the
    total number of unique players who have achieved this feat.
    """
    # Data source: MLB, Baseball-Reference, Wikipedia
    # Each dictionary represents a Triple Crown season win.
    triple_crown_wins = [
        # National League
        {'year': 1878, 'player': 'Paul Hines', 'league': 'NL'},
        {'year': 1894, 'player': 'Hugh Duffy', 'league': 'NL'},
        {'year': 1922, 'player': 'Rogers Hornsby', 'league': 'NL'},
        {'year': 1925, 'player': 'Rogers Hornsby', 'league': 'NL'},
        {'year': 1933, 'player': 'Chuck Klein', 'league': 'NL'},
        {'year': 1937, 'player': 'Joe Medwick', 'league': 'NL'},
        # American League
        {'year': 1901, 'player': 'Nap Lajoie', 'league': 'AL'},
        {'year': 1909, 'player': 'Ty Cobb', 'league': 'AL'},
        {'year': 1933, 'player': 'Jimmie Foxx', 'league': 'AL'},
        {'year': 1934, 'player': 'Lou Gehrig', 'league': 'AL'},
        {'year': 1937, 'player': 'Joe DiMaggio', 'league': 'AL'},
        {'year': 1942, 'player': 'Ted Williams', 'league': 'AL'},
        {'year': 1947, 'player': 'Ted Williams', 'league': 'AL'},
        {'year': 1956, 'player': 'Mickey Mantle', 'league': 'AL'},
        {'year': 1966, 'player': 'Frank Robinson', 'league': 'AL'},
        {'year': 1967, 'player': 'Carl Yastrzemski', 'league': 'AL'},
        {'year': 2012, 'player': 'Miguel Cabrera', 'league': 'AL'},
    ]

    # Separate winners by league
    al_winners = [win for win in triple_crown_wins if win['league'] == 'AL']
    nl_winners = [win for win in triple_crown_wins if win['league'] == 'NL']

    # Get the count of wins for each league
    num_al_wins = len(al_winners)
    num_nl_wins = len(nl_winners)
    
    # Calculate total instances
    total_instances = num_al_wins + num_nl_wins

    # Get the set of unique player names to find the number of unique winners
    unique_players = set(win['player'] for win in triple_crown_wins)
    num_unique_winners = len(unique_players)

    print(f"Number of Triple Crown wins in the American League: {num_al_wins}")
    print(f"Number of Triple Crown wins in the National League: {num_nl_wins}")
    print(f"Total Triple Crown instances in MLB history: {num_al_wins} + {num_nl_wins} = {total_instances}")
    print(f"\nThe number of unique players to win the Triple Crown is: {num_unique_winners}")

if __name__ == "__main__":
    count_triple_crown_winners()