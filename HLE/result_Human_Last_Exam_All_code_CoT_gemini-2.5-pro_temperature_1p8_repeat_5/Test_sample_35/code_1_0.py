import textwrap

def poker_bubble_analysis():
    """
    Analyzes a poker hand scenario on the money bubble and determines the correct action.
    """
    
    # --- 1. Define the Scenario ---
    position = "UTG+1 (Under the Gun + 1)"
    stack_size_bb = 16
    tournament_stage = "Near the money bubble"
    
    # --- 2. Define the Hand Choices ---
    hand_choices = {
        'A': 'QJs',
        'B': 'None of these',
        'C': '99',
        'D': 'AJo',
        'E': 'AKo'
    }

    # --- 3. Print the Analysis ---
    print(f"--- Poker Scenario Analysis ---")
    print(f"Position: {position}")
    print(f"Stack Size: {stack_size_bb} big blinds")
    print(f"Stage: {tournament_stage}\n")

    print("--- Evaluating the Hand Options ---\n")

    # Analysis for QJs
    analysis_qjs = """
    A. QJs (Queen-Jack suited): This is a playable hand in some situations, but shoving it from an early position like UTG+1 is very risky. With many players left to act, the probability of being called by a dominating hand (like AQ, KQ, AK, TT+) is high. The intense ICM pressure of the money bubble makes this a clear fold. You want to avoid marginal confrontations.
    """
    print(textwrap.dedent(analysis_qjs))
    
    # Analysis for 99
    analysis_99 = """
    C. 99 (Pocket Nines): This is a strong pocket pair. It is a standard all-in for 16bb in most tournament situations. It has good equity against the kinds of hands that might call (it's a 'coin flip' against two overcards like AK). While it's a very strong candidate for a shove, some extremely tight bubble situations could make it a borderline decision. It is, however, a very reasonable shove.
    """
    print(textwrap.dedent(analysis_99))

    # Analysis for AJo
    analysis_ajo = """
    D. AJo (Ace-Jack offsuit): This hand is strong, but it's also easily dominated by hands like AQ, AK, and JJ+, which are a significant part of an opponent's calling range. Shoving AJo from early position on the bubble exposes you to a high risk of being called and crushed. This makes it a likely fold in this specific spot.
    """
    print(textwrap.dedent(analysis_ajo))
    
    # Analysis for AKo
    analysis_ako = """
    E. AKo (Ace-King offsuit): This is a premium hand. It is one of the strongest starting hands in No-Limit Hold'em. It dominates a vast portion of your opponents' potential hands. You should always be happy to get all your chips in the middle with AK pre-flop for a 16bb stack. This is the most non-debatable, standard, and profitable shove among all the choices, regardless of bubble pressure.
    """
    print(textwrap.dedent(analysis_ako))

    # --- 4. Conclusion ---
    conclusion = """
    --- Conclusion ---
    Given the options, multiple hands could be considered for a shove in a vacuum (99, AJo, AKo). However, the combination of an early position (UTG+1) and severe bubble pressure requires a very tight and strong shoving range.
    
    - QJs and AJo are too weak and prone to domination.
    - 99 is a strong hand and a profitable shove.
    - AKo is an elite, premium hand.
    
    While 99 is a good shove, AKo is the strongest and most unambiguously correct hand to jam. It sits at the top of any shoving range and is a profitable all-in under any circumstances in this scenario.
    
    The best hand to jam is AKo.
    """
    print(textwrap.dedent(conclusion))

poker_bubble_analysis()
<<<E>>>