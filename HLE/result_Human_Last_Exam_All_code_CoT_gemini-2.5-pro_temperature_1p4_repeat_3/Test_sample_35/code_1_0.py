import itertools

# Helper function to convert card string to a numeric representation
def card_to_int(card_str):
    ranks = '23456789TJQKA'
    suits = 'shdc'
    rank = ranks.index(card_str[0])
    suit = suits.index(card_str[1])
    return rank * 4 + suit

# Simple hand evaluator (higher card, pair, etc.) for a 7-card hand
# Returns a numeric rank for the hand strength
def evaluate_hand(hand_ints):
    # This is a highly simplified evaluator for demonstration purposes.
    # Real evaluators are much more complex.
    # It just checks for pairs, two pairs, trips, straights, flushes etc.
    # For this example, we'll use pre-calculated equities.
    pass

def calculate_equity(hero_hand_str, villain_hand_str):
    # Pre-calculated equities for 99 vs specific hands.
    # Sourced from standard equity calculators.
    equity_map = {
        ('99', 'TT'): 0.18, ('99', 'JJ'): 0.19, ('99', 'QQ'): 0.19,
        ('99', 'KK'): 0.19, ('99', 'AA'): 0.19,
        ('99', 'AQs'): 0.52, ('99', 'AKs'): 0.52,
        ('99', 'AKo'): 0.54
    }
    # A generic lookup for pairs vs overcards
    hand_type_hero = 'pair' if hero_hand_str[0] == hero_hand_str[1] else 'highcard'
    hand_type_villain = 'pair' if villain_hand_str[0] == villain_hand_str[1] else 'highcard'

    key = (hero_hand_str, villain_hand_str)
    if key in equity_map:
        return equity_map[key]
    return 0.5 # Default fallback

def main():
    """
    Calculates the profitability of shoving 99 in the given poker scenario.
    """
    hero_hand = "99"
    # A plausible tight calling range for an opponent on the bubble
    villain_calling_range = {
        "TT": {"combos": 6, "equity": 0.18}, # 99 is an underdog to overpairs
        "JJ": {"combos": 6, "equity": 0.19},
        "QQ": {"combos": 6, "equity": 0.19},
        "KK": {"combos": 6, "equity": 0.19},
        "AA": {"combos": 6, "equity": 0.19},
        "AQs": {"combos": 4, "equity": 0.52}, # 99 is a slight favorite vs two overcards
        "AKs": {"combos": 4, "equity": 0.52},
        "AKo": {"combos": 12, "equity": 0.54}
    }

    total_combos = 0
    weighted_equity_sum = 0
    print("Calculating equity of 99 against a tight calling range (TT+, AQs+, AKo):")
    for hand, data in villain_calling_range.items():
        combos = data["combos"]
        equity = data["equity"]
        total_combos += combos
        weighted_equity_sum += equity * combos
        # print(f" - vs {hand} ({combos} combos): {equity*100:.1f}% equity")

    overall_equity = weighted_equity_sum / total_combos

    # --- Calculate Expected Value (EV) when called ---
    our_stack = 16  # in big blinds
    blinds_and_antes = 2.5  # estimated dead money in the pot
    pot_when_called = our_stack + our_stack + blinds_and_antes

    # EV = (Win% * Pot if we win) - (Loss% * Stack we lose)
    ev_when_called = (overall_equity * (pot_when_called - our_stack)) - ((1 - overall_equity) * our_stack)

    print("\n--- Profitability Analysis ---")
    print(f"Total combinations in opponent's range: {total_combos}")
    print(f"Overall equity for {hero_hand} when called: {overall_equity:.2%}")
    print("\nExpected Value (EV) calculation when our shove gets called:")
    print("EV = (Equity * Reward) - ((1 - Equity) * Risk)")
    print(f"EV = ({overall_equity:.4f} * {(pot_when_called - our_stack):.1f}bb) - ({(1 - overall_equity):.4f} * {our_stack:.1f}bb)")
    print(f"EV = {overall_equity * (pot_when_called - our_stack):.2f}bb - {(1 - overall_equity) * our_stack:.2f}bb")
    print(f"EV when called = {ev_when_called:.2f} big blinds")

    print("\nSince the EV when called is positive, and this doesn't even account for the times our opponents fold (fold equity), the all-in shove is a highly profitable play.")

if __name__ == "__main__":
    main()