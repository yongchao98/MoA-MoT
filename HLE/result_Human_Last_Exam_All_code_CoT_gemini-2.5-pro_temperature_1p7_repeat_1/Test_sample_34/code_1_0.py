def find_worst_suited_jack_to_open():
    """
    Determines the worst suited jack to open from the button in a 100BB rake-free game.

    Poker theory dictates that the button is the strongest pre-flop position,
    allowing for a very wide opening range. Game Theory Optimal (GTO) charts for
    a 100BB deep game indicate that all suited jacks are profitable opens from this position,
    especially in a rake-free environment which makes marginal hands more playable.

    The suited jacks are ranked by the value of their second card (kicker).
    The hand 'J2s' (Jack-Two suited) has the lowest possible kicker.
    Since all suited jacks are considered standard opens in this scenario,
    the "worst" one is the one with the lowest rank.
    """
    
    # The hand is abbreviated as the higher card, the lower card, and 's' for suited.
    worst_suited_jack = "J2s"
    
    print("Based on Game Theory Optimal (GTO) strategy for a 100BB rake-free cash game:")
    print("The worst suited Jack that you should open-raise from the button position is J2s.")
    print("The hand is represented by the characters:")
    
    # Per instructions, outputting each character in the final "equation" (in this case, the hand name)
    print(f"Character 1: {worst_suited_jack[0]}")
    print(f"Character 2: {worst_suited_jack[1]}")
    print(f"Character 3: {worst_suited_jack[2]}")

find_worst_suited_jack_to_open()