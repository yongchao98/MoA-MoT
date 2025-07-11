def count_triple_crown_winners():
    """
    Calculates and displays the number of unique MLB Triple Crown winners.
    """
    # This list represents all recognized batting Triple Crown winning seasons in
    # MLB history, including the American Association (AA), which was a major
    # league in the 19th century.
    triple_crown_seasons = [
        {'player': 'Paul Hines', 'year': 1878, 'league': 'NL'},
        {'player': 'Tip O\'Neill', 'year': 1887, 'league': 'AA'},
        {'player': 'Hugh Duffy', 'year': 1894, 'league': 'NL'},
        {'player': 'Nap Lajoie', 'year': 1901, 'league': 'AL'},
        {'player': 'Ty Cobb', 'year': 1909, 'league': 'AL'},
        {'player': 'Heinie Zimmerman', 'year': 1912, 'league': 'NL'},
        {'player': 'Rogers Hornsby', 'year': 1922, 'league': 'NL'},
        {'player': 'Rogers Hornsby', 'year': 1925, 'league': 'NL'},
        {'player': 'Jimmie Foxx', 'year': 1933, 'league': 'AL'},
        {'player': 'Chuck Klein', 'year': 1933, 'league': 'NL'},
        {'player': 'Lou Gehrig', 'year': 1934, 'league': 'AL'},
        {'player': 'Joe Medwick', 'year': 1937, 'league': 'NL'},
        {'player': 'Ted Williams', 'year': 1942, 'league': 'AL'},
        {'player': 'Ted Williams', 'year': 1947, 'league': 'AL'},
        {'player': 'Mickey Mantle', 'year': 1956, 'league': 'AL'},
        {'player': 'Frank Robinson', 'year': 1966, 'league': 'AL'},
        {'player': 'Carl Yastrzemski', 'year': 1967, 'league': 'AL'},
        {'player': 'Miguel Cabrera', 'year': 2012, 'league': 'AL'}
    ]

    # Use a set to get a list of unique players, then sort it alphabetically.
    # The question asks for the number of "winners", which implies unique people.
    unique_winners = sorted(list(set(season['player'] for season in triple_crown_seasons)))
    num_winners = len(unique_winners)

    print("To find the total number of unique MLB Triple Crown winners, we add each player who has won the award:")
    
    equation_parts = []
    # Display each unique player as '1' in our count.
    for winner in unique_winners:
        # Each player counts as 1 winner, regardless of how many times they won.
        print(f"1 ({winner})")
        equation_parts.append('1')
        
    equation_str = " + ".join(equation_parts)

    print("\nFinal Equation:")
    print(f"{equation_str} = {num_winners}")

# Run the function to display the result.
count_triple_crown_winners()