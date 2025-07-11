def find_worst_suited_jack_to_open():
    """
    This function determines the worst suited Jack hand that should be open-raised
    from the button in a 100BB deep, rake-free cash game based on standard
    Game Theory Optimal (GTO) poker strategy.

    In GTO, all suited Jacks (from AJs down to J2s) are standard opens from the
    button. Therefore, the "worst" one is the one with the lowest kicker.
    """
    
    # The components of the hand
    high_card = "J"
    low_card_rank = 2  # The lowest possible kicker for a Jack
    suit_indicator = "s"  # 's' for suited
    
    # Construct the abbreviated hand name
    # The number '2' is part of the final hand name
    worst_hand = f"{high_card}{low_card_rank}{suit_indicator}"
    
    print("Based on standard GTO ranges, the worst suited jack to open from the button is:")
    print(worst_hand)

if __name__ == "__main__":
    find_worst_suited_jack_to_open()