def solve_poker_question():
    """
    This script determines and prints the worst suited Jack a player should
    open-raise from the button in a 100BB rake-free cash game.
    """

    # In poker, a hand is often abbreviated. 'J' stands for Jack. A number
    # represents the other card's rank. 's' means the cards are 'suited'
    # (i.e., of the same suit, like two hearts).

    # Game Theory Optimal (GTO) charts show the profitable hands to play.
    # From the button position at 100BB depth, the opening range is very wide.
    # This range includes all suited Jack hands, from the strongest (AJs)
    # to the weakest.

    # The "worst" hand to open is the one in this range with the lowest
    # expected value. We must identify the weakest suited Jack.

    # J2s is the weakest because it has the lowest kicker (the 2) and the
    # worst connectivity for making straights.
    
    # We will represent the hand by its components:
    # High card: 'J'
    # Low card (number): 2
    # Suit status: 's' for suited

    high_card = "J"
    low_card_number = 2
    suit_status = "s"

    # Combine the parts to form the standard hand abbreviation.
    worst_suited_jack_hand = f"{high_card}{low_card_number}{suit_status}"

    print("In a 100BB rake-free cash game, the GTO opening range from the button is very wide.")
    print("This range includes all suited Jack hands because of their flush potential and positional advantage.")
    print("The worst of these hands, and therefore the marginal open, is the one with the lowest kicker and weakest connectivity.")
    print("-" * 30)
    print("The worst suited Jack you should open from the button is:")
    print(worst_suited_jack_hand)

solve_poker_question()