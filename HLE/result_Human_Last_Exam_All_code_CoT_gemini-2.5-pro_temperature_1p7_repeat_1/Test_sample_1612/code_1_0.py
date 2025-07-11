def count_mlb_triple_crowns():
    """
    Calculates and prints the total number of MLB Triple Crown instances (batting and pitching).
    The data includes winners from the modern era (post-1900).
    """

    # Data for Batting Triple Crown winners
    batting_winners = [
        {"player": "Nap Lajoie", "year": 1901},
        {"player": "Ty Cobb", "year": 1909},
        {"player": "Heinie Zimmerman", "year": 1912},
        {"player": "Rogers Hornsby", "year": 1922},
        {"player": "Rogers Hornsby", "year": 1925},
        {"player": "Jimmie Foxx", "year": 1933},
        {"player": "Chuck Klein", "year": 1933},
        {"player": "Lou Gehrig", "year": 1934},
        {"player": "Joe Medwick", "year": 1937},
        {"player": "Ted Williams", "year": 1942},
        {"player": "Ted Williams", "year": 1947},
        {"player": "Mickey Mantle", "year": 1956},
        {"player": "Frank Robinson", "year": 1966},
        {"player": "Carl Yastrzemski", "year": 1967},
        {"player": "Miguel Cabrera", "year": 2012},
    ]

    # Data for Pitching Triple Crown winners
    pitching_winners = [
        {"player": "Cy Young", "year": 1901},
        {"player": "Rube Waddell", "year": 1905},
        {"player": "Christy Mathewson", "year": 1905},
        {"player": "Christy Mathewson", "year": 1908},
        {"player": "Walter Johnson", "year": 1913},
        {"player": "Pete Alexander", "year": 1915},
        {"player": "Pete Alexander", "year": 1916},
        {"player": "Walter Johnson", "year": 1918},
        {"player": "Hippo Vaughn", "year": 1918},
        {"player": "Pete Alexander", "year": 1920},
        {"player": "Walter Johnson", "year": 1924},
        {"player": "Dazzy Vance", "year": 1924},
        {"player": "Lefty Grove", "year": 1930},
        {"player": "Lefty Grove", "year": 1931},
        {"player": "Lefty Gomez", "year": 1934},
        {"player": "Lefty Gomez", "year": 1937},
        {"player": "Bucky Walters", "year": 1939},
        {"player": "Bob Feller", "year": 1940},
        {"player": "Hal Newhouser", "year": 1945},
        {"player": "Sandy Koufax", "year": 1963},
        {"player": "Sandy Koufax", "year": 1965},
        {"player": "Sandy Koufax", "year": 1966},
        {"player": "Steve Carlton", "year": 1972},
        {"player": "Dwight Gooden", "year": 1985},
        {"player": "Roger Clemens", "year": 1997},
        {"player": "Roger Clemens", "year": 1998},
        {"player": "Pedro Martinez", "year": 1999},
        {"player": "Randy Johnson", "year": 2002},
        {"player": "Johan Santana", "year": 2006},
        {"player": "Jake Peavy", "year": 2007},
        {"player": "Justin Verlander", "year": 2011},
        {"player": "Clayton Kershaw", "year": 2011},
        {"player": "Shane Bieber", "year": 2020},
    ]

    num_batting_crowns = len(batting_winners)
    num_pitching_crowns = len(pitching_winners)
    total_crowns = num_batting_crowns + num_pitching_crowns

    print(f"To find the total number of MLB Triple Crowns, we sum the batting and pitching instances.")
    print(f"Total Triple Crowns = (Batting Triple Crowns) + (Pitching Triple Crowns)")
    print(f"Total = {num_batting_crowns} + {num_pitching_crowns} = {total_crowns}")

if __name__ == "__main__":
    count_mlb_triple_crowns()