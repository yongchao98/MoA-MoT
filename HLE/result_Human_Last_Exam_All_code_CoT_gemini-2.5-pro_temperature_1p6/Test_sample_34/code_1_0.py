def find_worst_suited_jack():
    """
    Identifies the worst suited Jack to open from the button in a standard 100BB cash game.
    
    In Game Theory Optimal (GTO) play from the button, all suited Jacks are considered
    profitable opens due to the significant positional advantage. The "worst" hand in
    this category is the one with the lowest kicker, as it has the least equity
    and is most easily dominated.
    """
    
    # List of suited Jack hands, ranked from strongest to weakest
    suited_jacks = ["JTs", "J9s", "J8s", "J7s", "J6s", "J5s", "J4s", "J3s", "J2s"]
    
    # The last hand in the ranked list is the weakest one that is still a standard open
    worst_suited_jack = suited_jacks[-1]
    
    print(f"The worst suited jack you should open from the button is: {worst_suited_jack}")

find_worst_suited_jack()