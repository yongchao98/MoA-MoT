import math

def solve_poker_strategy():
    """
    Calculates the optimal GTO strategy for the given poker scenario.
    """
    # 1. Define the Game State
    pot = 10
    stack = 1000
    # Our range is AA (1 value combo) and QQ-33 (10 bluff candidate combos).
    # Villain's range is KK (a pure bluff-catcher).

    # 2. The GTO Approach
    # In a GTO strategy, we must make our opponent indifferent to certain actions,
    # which in turn makes us unexploitable.
    # We have two types of hands: value (AA) and bluffs (QQ-33).
    # The core of the strategy is to bet with our value hand and balance it
    # by bluffing with a specific frequency of our weaker hands.

    # 3. Making our Bluffs Indifferent
    # Our bluffs (QQ-33) must be indifferent between betting and checking.
    # EV(Check) = 0. If we check, the villain knows we are weak, will bet, and we will fold.
    # EV(Bluff Bet) must also be 0.
    # EV(Bluff Bet) = (Prob_Villain_Folds * Pot) - (Prob_Villain_Calls * Bet_Size)
    # Let S be the Bet_Size and p_call be the villain's calling frequency.
    # (1 - p_call) * pot - p_call * S = 0
    # This means the villain must call with a frequency p_call = pot / (pot + S).
    # This is the Minimum Defense Frequency (MDF).

    # 4. Maximizing our Value Hand's EV
    # Since our bluffs are constructed to have an EV of 0, our total EV is determined
    # by the EV of our value hand, AA. We should choose a bet size S that maximizes this.
    # EV(AA Bet) = (Prob_Villain_Calls * Total_Pot_Won) + (Prob_Villain_Folds * Pot_Won)
    # When villain calls, we win pot + S. When he folds, we win the pot.
    # EV(AA Bet) = (p_call * (pot + S)) + ((1 - p_call) * pot)

    # We substitute p_call = pot / (pot + S):
    # EV(AA Bet) = (pot / (pot + S) * (pot + S)) + ((1 - pot / (pot + S)) * pot)
    # EV(AA Bet) = pot + (S / (pot + S) * pot)

    # 5. Determine the Optimal Bet Size
    # To maximize EV(AA Bet), we need to maximize the term: pot + pot * S / (pot + S).
    # Since 'pot' is constant, we just need to maximize S / (pot + S).
    # This function increases as S increases. Therefore, to maximize our EV,
    # we should bet the largest amount possible.
    optimal_bet_size = stack

    # 6. Formulate the Final Strategy for AA
    # The optimal strategy for our value hand, AA, is to always bet the maximum amount.
    action = "BET"
    sizing = optimal_bet_size
    frequency = 100 # We always take the highest EV line with our value hands.

    # Round sizing and frequency to the nearest even number as per instructions.
    # 1000 is already even.
    if sizing % 2 != 0:
        sizing = int(round(sizing / 2.0) * 2.0)
    # 100 is already even.
    if frequency % 2 != 0:
        frequency = int(round(frequency / 2.0) * 2.0)

    # Print the final strategy for holding AA.
    # The final equation is the representation of the strategy itself.
    print(f"{action} {sizing} {frequency}%")


solve_poker_strategy()
<<<BET 1000 100%>>>