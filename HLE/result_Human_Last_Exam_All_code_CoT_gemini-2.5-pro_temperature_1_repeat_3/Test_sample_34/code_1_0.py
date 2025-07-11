import collections

def find_worst_suited_jack():
    """
    Determines the worst suited Jack to open from the button in a 100BB rake-free cash game.

    In poker, the value of a suited hand is strongly influenced by its kicker.
    A standard GTO (Game Theory Optimal) opening range from the button is very wide,
    especially in a rake-free environment. This range includes all suited Jack hands.

    The "worst" hand to open is the one with the lowest profitability, which is typically
    the one with the lowest kicker. A lower kicker means weaker flush potential (a J-high flush
    is more likely to be beaten) and weaker straight potential.
    """

    # Define the card ranks in descending order of strength (Ace is high)
    # The kicker is the second card in a JXs hand.
    ranks = ['A', 'K', 'Q', 'T', '9', '8', '7', '6', '5', '4', '3', '2']

    # Create a list of all suited Jack hands
    suited_jacks = [f"J{rank}s" for rank in ranks]

    print("List of all suited Jack hands, from best to worst:")
    print(suited_jacks)
    print("-" * 30)

    # The worst hand is the one with the lowest kicker, which is the last in the ordered list.
    worst_hand = suited_jacks[-1]

    # Deconstruct the hand name as requested by the prompt to show the "equation"
    card1 = worst_hand[0]
    card2 = worst_hand[1]
    hand_type = worst_hand[2]

    print(f"The 'equation' for the hand is combining its components:")
    print(f"High Card: {card1}")
    print(f"Kicker: {card2}")
    print(f"Type: {hand_type} (suited)")
    print("-" * 30)
    print(f"Based on poker theory, the worst suited Jack that is still a standard open from the button is: {worst_hand}")

find_worst_suited_jack()