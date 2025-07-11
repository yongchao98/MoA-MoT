import sys

def analyze_shove_decision():
    """
    Analyzes the correct hand to shove UTG+1 with 16bb on the bubble.
    """

    # --- Step 1: Define the Scenario & Shoving Range ---
    position = "UTG+1"
    stack_bb = 16
    situation = "Near the money bubble"

    # On the bubble from early position, we must tighten our shoving range
    # compared to a standard push/fold chart. We want to reduce risk and
    # shove hands that dominate potential calling ranges.
    # A strong, conservative range is: TT+, AQs+, AKo.
    # - Pairs: TT, JJ, QQ, KK, AA
    # - Suited Aces: AQs, AKs
    # - Offsuit Hands: AKo
    # Let's represent the lowest hand of each category required for a shove.
    min_pair_rank_index = "23456789TJQKA".index('T')
    min_suited_ace_rank_index = "23456789TJQKA".index('Q')
    
    # --- Step 2: Define Hand Options ---
    options = {
        "A": "QJs",
        "C": "99",
        "D": "AJo",
        "E": "AKo"
    }
    
    print(f"Scenario: {stack_bb}bb stack in {position}, {situation}.")
    print("Analysis based on a conservative shoving range of TT+, AQs+, AKo.")
    print("-" * 60)

    # --- Step 3: Evaluate each hand ---
    # Hand A: QJs
    hand_a = options["A"]
    is_shove_a = False
    print(f"Hand A ({hand_a}): Is it a shove?")
    print(f"Evaluation: {is_shove_a}. QJs is too weak to shove from this position on the bubble.")
    print("-" * 60)

    # Hand C: 99
    hand_c = options["C"]
    hand_c_rank_index = "23456789TJQKA".index('9')
    # The "equation" here is a comparison of rank indices.
    # Is the rank of 9 (index 7) >= the rank of T (index 8)?
    is_shove_c = hand_c_rank_index >= min_pair_rank_index
    print(f"Hand C ({hand_c}): Is it a shove?")
    print(f"Evaluation: Rank of '9' ({hand_c_rank_index}) >= Minimum Pair Rank 'T' ({min_pair_rank_index})? -> {is_shove_c}")
    print("While 99 is a good hand, it falls just outside a conservative bubble range of TT+.")
    print("-" * 60)

    # Hand D: AJo
    hand_d = options["D"]
    is_shove_d = False
    print(f"Hand D ({hand_d}): Is it a shove?")
    print(f"Evaluation: {is_shove_d}. AJo is dominated by many hands that will call (AK, AQ, TT+). It's not in our AQs+ or AKo range.")
    print("-" * 60)

    # Hand E: AKo
    hand_e = options["E"]
    # For AKo, the check is simple: is the hand literally "AKo"?
    is_shove_e = (hand_e == "AKo")
    print(f"Hand E ({hand_e}): Is it a shove?")
    print(f"Evaluation: Is the hand 'AKo'? -> {is_shove_e}")
    print("AKo is a premium hand that dominates calling ranges. It is a clear and standard shove.")
    print("-" * 60)

    print("\nConclusion: AKo is the strongest and most appropriate hand to jam in this specific scenario.")

if __name__ == '__main__':
    analyze_shove_decision()