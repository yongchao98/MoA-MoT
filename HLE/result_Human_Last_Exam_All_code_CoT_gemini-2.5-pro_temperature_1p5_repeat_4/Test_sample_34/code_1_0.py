def find_worst_suited_jack_open():
    """
    Identifies the worst suited jack to open-raise from the button in a
    standard 100BB NLHE cash game based on GTO principles.

    In poker theory, the button is the most profitable position, allowing for a
    very wide opening range. Game Theory Optimal (GTO) charts show that all
    suited jacks (from JAs down to J2s) are standard opens from this position.

    Therefore, the "worst" hand in this category that is still a profitable
    open-raise is the one with the lowest kicker.
    """
    
    # List of suited jacks from best to worst
    suited_jacks = ["JAs", "JKs", "JQs", "JTs", "J9s", "J8s", "J7s", "J6s", "J5s", "J4s", "J3s", "J2s"]
    
    # According to GTO charts, J2s is the last and therefore "worst" suited jack
    # that is included in a standard button opening range.
    worst_suited_jack_open = suited_jacks[-1]
    
    print(f"Based on standard GTO poker opening ranges, the worst suited jack you should open from the button is: {worst_suited_jack_open}")

find_worst_suited_jack_open()