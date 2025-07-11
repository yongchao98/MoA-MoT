def find_worst_opening_hand():
    """
    This function determines the worst suited Jack to open-raise from the button
    in a standard 100BB cash game, based on established poker theory.
    """

    # In poker hand notation, a hand is defined by its cards.
    # We are looking for a hand containing a Jack.
    main_card = "J"

    # A standard button opening range includes all suited Jacks.
    # The strength of these hands is determined by the second card, or "kicker".
    # The suited Jacks in this range are JTs, J9s, J8s, J7s, J6s, J5s, J4s, J3s, and J2s.
    # The "worst" of these hands is the one with the lowest kicker.
    worst_kicker = "2"

    # The 's' indicates that the two cards are of the same suit.
    suit_type = "s"

    # We combine these components to form the final hand abbreviation.
    # The final equation is effectively: main_card + worst_kicker + suit_type
    final_hand = main_card + worst_kicker + suit_type

    print(f"The main card is: {main_card}")
    print(f"The lowest kicker in the standard opening range for this hand type is: {worst_kicker}")
    print(f"The modifier for a suited hand is: {suit_type}")
    print("---")
    print(f"The worst suited jack you should open from the button is: {final_hand}")

find_worst_opening_hand()