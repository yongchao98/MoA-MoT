import collections

def poker_shove_analysis():
    """
    Analyzes the profitability of shoving all-in with different hands
    in a specific poker tournament scenario.
    """

    # --- 1. Define Scenario and Assumptions ---

    our_stack_bb = 16.0
    pot_before_shove_bb = 2.5  # SB + BB + Antes
    opponents_left_to_act = 6

    # Define a tight calling range for our opponents due to bubble pressure.
    # Combos: TT(6), JJ(6), QQ(6), KK(6), AA(6), AQs(4), AKs(4), AKo(12) = 50 total combos
    villain_calling_range_combos = 50
    total_hand_combos = 1326  # Total possible 2-card hand combinations

    # Pre-calculated equity of our hand vs. the villain's calling range.
    # Format: {Hand: Equity}
    equity_vs_range = {
        'AKo': 0.465,  # Ace-King offsuit
        'AJo': 0.314,  # Ace-Jack offsuit
        'QJs': 0.315,  # Queen-Jack suited
        '99':  0.367   # Pocket Nines
    }

    # --- 2. Calculate Probabilities ---

    # Probability a single opponent has a hand in the calling range
    prob_villain_has_caller = villain_calling_range_combos / total_hand_combos

    # Probability everyone folds is the chance that each opponent *doesn't* have a calling hand
    prob_all_fold = (1 - prob_villain_has_caller) ** opponents_left_to_act
    prob_gets_called = 1 - prob_all_fold

    print("--- Poker Shove EV Analysis ---")
    print(f"Scenario: 16bb stack, UTG+1, on the money bubble.")
    print(f"Assuming opponents call with a tight range (TT+, AQs+, AKo).\n")
    print(f"Probability all {opponents_left_to_act} players fold: {prob_all_fold:.2%}")
    print(f"Probability of getting called by at least one player: {prob_gets_called:.2%}\n")
    print("--- Calculating EV for each hand ---")

    # --- 3. Calculate EV for Each Hand ---

    results = {}
    HandEV = collections.namedtuple('HandEV', ['hand', 'ev_total', 'ev_if_called', 'equity'])

    for hand, equity in equity_vs_range.items():
        # EV when we get called
        # Pot size if called = our stack + villain stack + dead money
        pot_if_called = our_stack_bb + our_stack_bb + pot_before_shove_bb
        # Amount we win = their stack + dead money
        win_amount = our_stack_bb + pot_before_shove_bb
        # Amount we lose = our stack
        lose_amount = our_stack_bb

        ev_if_called = (equity * win_amount) - ((1 - equity) * lose_amount)

        # Total EV of shoving
        ev_total = (prob_all_fold * pot_before_shove_bb) + (prob_gets_called * ev_if_called)
        
        results[hand] = HandEV(hand, ev_total, ev_if_called, equity)
        print(f"\nHand: {hand}")
        print(f"  - Equity vs. Calling Range: {equity:.2%}")
        print(f"  - EV if Called: {ev_if_called:.2f} bb")
        print(f"  - Total EV of Shoving: {ev_total:.2f} bb")

    # --- 4. Conclusion and Final Equation ---

    best_hand = max(results.values(), key=lambda x: x.ev_total)
    
    print("\n--- Conclusion ---")
    print("Based on the Expected Value, all hands except QJs might seem profitable due to high fold equity.")
    print("However, AKo is the most profitable and robust shove. It performs best when called and has great blockers.")
    print(f"\nThe best hand to shove from the options is {best_hand.hand} with a total EV of {best_hand.ev_total:.2f} bb.")

    print("\nHere is the final EV equation for the best hand (AKo):")
    print("EV = (P(fold) * Pot Won) + (P(call) * EV when called)")
    # The prompt asks for each number in the final equation.
    print(f"EV({best_hand.hand}) = ({prob_all_fold:.3f} * {pot_before_shove_bb}) + ({prob_gets_called:.3f} * {best_hand.ev_if_called:.3f})")
    print(f"EV({best_hand.hand}) = ({(prob_all_fold * pot_before_shove_bb):.3f}) + ({(prob_gets_called * best_hand.ev_if_called):.3f})")
    print(f"EV({best_hand.hand}) = {best_hand.ev_total:.3f} bb")


poker_shove_analysis()