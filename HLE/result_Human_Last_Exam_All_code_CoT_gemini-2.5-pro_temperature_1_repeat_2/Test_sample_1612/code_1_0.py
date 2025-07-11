def count_mlb_triple_crown_winners():
    """
    This function calculates and prints the number of MLB Triple Crown winners.
    It counts unique players for both batting and pitching Triple Crowns and then sums them.
    """
    # Data for Batting Triple Crown winners (Name, Year).
    # Some players like Rogers Hornsby and Ted Williams won it multiple times.
    batting_winners_data = [
        ("Paul Hines", 1878), ("Tip O'Neill", 1887), ("Hugh Duffy", 1894),
        ("Nap Lajoie", 1901), ("Ty Cobb", 1909), ("Rogers Hornsby", 1922),
        ("Rogers Hornsby", 1925), ("Chuck Klein", 1933), ("Jimmie Foxx", 1933),
        ("Lou Gehrig", 1934), ("Joe Medwick", 1937), ("Ted Williams", 1942),
        ("Ted Williams", 1947), ("Mickey Mantle", 1956), ("Frank Robinson", 1966),
        ("Carl Yastrzemski", 1967), ("Miguel Cabrera", 2012)
    ]

    # Data for Pitching Triple Crown winners (Name, Year).
    # Many pitchers have won multiple times (e.g., Walter Johnson, Sandy Koufax).
    pitching_winners_data = [
        ("Tommy Bond", 1877), ("Old Hoss Radbourn", 1884), ("Tim Keefe", 1888),
        ("John Clarkson", 1889), ("Amos Rusie", 1894), ("Cy Young", 1901),
        ("Rube Waddell", 1905), ("Christy Mathewson", 1905), ("Christy Mathewson", 1908),
        ("Walter Johnson", 1913), ("Walter Johnson", 1918), ("Walter Johnson", 1924),
        ("Pete Alexander", 1915), ("Pete Alexander", 1916), ("Pete Alexander", 1920),
        ("Dazzy Vance", 1924), ("Lefty Grove", 1930), ("Lefty Grove", 1931),
        ("Lefty Gomez", 1934), ("Lefty Gomez", 1937), ("Bucky Walters", 1939),
        ("Hal Newhouser", 1945), ("Sandy Koufax", 1963), ("Sandy Koufax", 1965),
        ("Sandy Koufax", 1966), ("Steve Carlton", 1972), ("Roger Clemens", 1997),
        ("Roger Clemens", 1998), ("Pedro Martínez", 1999), ("Randy Johnson", 2002),
        ("Johan Santana", 2006), ("Clayton Kershaw", 2011), ("Justin Verlander", 2011),
        ("Shane Bieber", 2020)
    ]

    # Use sets to find the number of unique players for each category.
    unique_batting_winners = set(player[0] for player in batting_winners_data)
    unique_pitching_winners = set(player[0] for player in pitching_winners_data)

    num_batting_winners = len(unique_batting_winners)
    num_pitching_winners = len(unique_pitching_winners)
    total_winners = num_batting_winners + num_pitching_winners

    print(f"Number of unique batting Triple Crown winners: {num_batting_winners}")
    print(f"Number of unique pitching Triple Crown winners: {num_pitching_winners}")
    print(f"Total number of MLB Triple Crown winners is the sum of both categories:")
    print(f"{num_batting_winners} + {num_pitching_winners} = {total_winners}")

count_mlb_triple_crown_winners()
<<<39>>>