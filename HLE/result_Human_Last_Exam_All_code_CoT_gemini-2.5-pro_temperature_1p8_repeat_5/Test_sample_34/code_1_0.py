def find_worst_suited_jack():
    """
    In a 100BB rake-free cash game, Game Theory Optimal (GTO) strategy advises
    opening a very wide range of hands from the button due to the positional advantage.
    
    This range includes all suited Jacks, which are ranked as follows from best to worst:
    JTs, J9s, J8s, J7s, J6s, J5s, J4s, J3s, J2s.
    
    The "worst" hand in this category is the one with the lowest kicker that is still
    considered a standard, profitable open.
    """
    
    # List of all suited Jacks that are a profitable open from the button.
    suited_jacks_range = ['JTs', 'J9s', 'J8s', 'J7s', 'J6s', 'J5s', 'J4s', 'J3s', 'J2s']
    
    # The last hand in the ranked list is the one with the lowest equity.
    worst_suited_jack = suited_jacks_range[-1]
    
    print(worst_suited_jack)

find_worst_suited_jack()