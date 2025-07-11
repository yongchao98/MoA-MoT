def find_worst_suited_jack_to_open():
    """
    This function determines the worst suited Jack to open from the button
    in a 100BB rake-free cash game based on GTO principles.
    """
    
    # List of all suited Jack hands.
    suited_jacks = ["JAs", "JKs", "JQs", "JTs", "J9s", "J8s", "J7s", "J6s", "J5s", "J4s", "J3s", "J2s"]
    
    print("Poker Hand Analysis: Worst Suited Jack to Open from the Button")
    print("="*60)
    print("Scenario: 100BB deep, Rake-Free, Open-raising from the Button (BTN).")
    print("\nReasoning:")
    print("1. Positional Advantage: The Button is the best position, allowing for the widest opening range.")
    print("2. Rake-Free Environment: With no rake, even marginal hands become profitable to open.")
    print("3. Evaluating Hands: We are looking for the suited Jack with the lowest Expected Value (EV) that is still a standard open.")
    print("\nComparison:")
    print("- Strong suited Jacks (JTs, JQs, JKs) are premium hands.")
    print("- Medium suited Jacks (J9s, J8s, J7s) have good flush and straight potential.")
    print("- Weak suited Jacks (J6s, J5s, J4s, J3s, J2s) are at the bottom of the range.")
    
    # Identify the worst hand based on poker theory.
    # The primary value is the flush draw.
    # The secondary value comes from straight draws and pairing the Jack.
    # J2s has the worst kicker and the least straight connectivity.
    # For example, J5s can help make a straight (A-2-3-4-5), while J2s has almost no meaningful straight potential.
    # Therefore, J2s has the lowest EV among all suited Jacks that are standard opens in this spot.
    
    worst_hand = "J2s"
    
    print(f"\nConclusion: The hand '{worst_hand}' has the weakest kicker and least connectivity, making it the worst suited Jack to open.")
    
    print("\nFinal Answer:")
    print(worst_hand)

find_worst_suited_jack_to_open()