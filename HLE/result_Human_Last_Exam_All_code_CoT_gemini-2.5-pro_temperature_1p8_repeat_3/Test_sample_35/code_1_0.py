import textwrap

def analyze_poker_hand():
    """
    Analyzes a poker hand scenario to determine the correct all-in shove.
    """
    # -- Scenario Details --
    position = "UTG+1 (Under The Gun + 1)"
    stack_size_bb = 16
    situation = "Near the money bubble"

    print("--- Poker Hand Analysis ---")
    print(f"Stack Size: {stack_size_bb} big blinds")
    print(f"Position: {position}")
    print(f"Situation: {situation}\n")

    analysis_header = "Analysis:"
    analysis_text = """
    This is a 'shove-or-fold' stack size from an early position. The crucial factor is the money bubble, which introduces significant ICM pressure. This means we must be much tighter with our shoving range than usual. We want to shove hands that have excellent equity and good 'blocker' effects (holding an Ace or King reduces the chances opponents have premium hands) when called.
    """
    print(analysis_header)
    print(textwrap.fill(analysis_text, 80))
    print("\n--- Evaluating Hand Options ---\n")

    # -- Hand Options Analysis --
    options = {
        'A': 'QJs (Queen-Jack suited)',
        'C': '99 (Pocket Nines)',
        'D': 'AJo (Ace-Jack offsuit)',
        'E': 'AKo (Ace-King offsuit)'
    }

    # Analysis for each hand
    print(f"A. Hand: {options['A']}")
    print("   Verdict: FOLD. This hand is too weak to shove from UTG+1 on a bubble. It is often dominated by the strong hands that will call you (e.g., AQ, KQ, JJ+).\n")

    print(f"D. Hand: {options['D']}")
    print("   Verdict: FOLD. While stronger than QJs, AJo is also a dangerous hand here. It's at the very bottom of a non-bubble shoving range and becomes a clear fold with ICM pressure. You will be dominated by AQ, AK when called.\n")
    
    print(f"C. Hand: {options['C']}")
    print("   Verdict: STANDARD SHOVE. Pocket nines is a strong hand and a profitable shove in this spot. Its main weakness is running into overpairs (TT, JJ, QQ+) or being a coin flip against hands like AK and AQ.\n")
    
    print(f"E. Hand: {options['E']}")
    print("   Verdict: BEST SHOVE. AKo is a premium hand. It dominates almost all other Ace-X and King-X hands. Its blocker effect is powerful, making it less likely you run into AA or KK. Even against a pocket pair like QQ, you have nearly 50% equity. It maintains its value extremely well even under heavy ICM pressure.\n")

    # -- Final Conclusion --
    print("--- Conclusion ---")
    conclusion_text = """
    Both 99 and AKo are considered standard shoves in this spot. However, AKo is the superior choice. It performs better against a tight bubble calling range (TT+, AQ+) and has significant blocker value. In extreme bubble situations where survival is paramount, some might even fold 99, but AKo almost always remains a profitable shove. Therefore, it is the best and most robust choice among the options.
    """
    print(textwrap.fill(conclusion_text, 80))

# Execute the analysis
if __name__ == "__main__":
    analyze_poker_hand()