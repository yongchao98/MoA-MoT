def solve_poker_strategy():
    """
    Calculates the optimal GTO poker strategy for the given river scenario,
    adhering to the constraint that the bet sizing and all action
    frequencies must be even integers.
    """
    pot = 10
    stack = 1000
    optimal_s = 0

    # We must find the largest even bet size 'S' that results in even
    # integer percentages for both hero's bluffs and villain's calls.
    # We loop downwards from the max possible bet to find the optimal S.
    for s in range(stack, 0, -2):
        # Villain's calling percentage is c_pct = 100 * Pot / (S + Pot)
        # This simplifies to 1000 / (S + 10).
        # For c_pct to be an integer, (S + 10) must be a divisor of 1000.
        # For c_pct to be even, the result of the division must be even.
        if 1000 % (s + pot) == 0:
            c_pct = 1000 / (s + pot)
            if c_pct % 2 == 0:
                # If villain's call % is even, hero's bluff % will also be
                # even, since Bluff % = 100% - Call %.
                # We have found the largest S that works.
                optimal_s = s
                break
    
    # Calculate final percentages based on the optimal bet size
    villain_call_pct = int(100 * pot / (optimal_s + pot))
    hero_bluff_pct = int(100 * optimal_s / (optimal_s + pot))
    hero_check_pct = 100 - hero_bluff_pct
    villain_fold_pct = 100 - villain_call_pct

    print("To find the optimal strategy, we solve for the largest bet size 'S' where all resulting strategic percentages are even numbers.")
    print(f"The optimal bet size found is S = ${optimal_s}.")
    print("\nThe calculations are as follows:")
    print("---")
    
    print("1. Villain's required calling frequency (to make Hero's QQ bluffs break-even):")
    print(f"   c = Pot / (S + Pot) = {pot} / ({optimal_s} + {pot}) = {villain_call_pct / 100.0:.2f}")
    print(f"   The percentage is {villain_call_pct}%, which is an even number.")
    print("")

    print("2. Hero's required bluffing frequency with QQ (to make Villain's KK calls break-even):")
    print(f"   y = S / (S + Pot) = {optimal_s} / ({optimal_s} + {pot}) = {hero_bluff_pct / 100.0:.2f}")
    print(f"   The percentage is {hero_bluff_pct}%, which is an even number.")
    print("---")

    # Format the final answer string
    final_answer_str = (
        f"A) BET ${optimal_s} 100% "
        f"B) BET ${optimal_s} {hero_bluff_pct}%, CHECK {hero_check_pct}% "
        f"C) CALL {villain_call_pct}%, FOLD {villain_fold_pct}%"
    )

    print("\nBased on this, the final optimal strategy is:")
    print(final_answer_str)

# Execute the function to print the solution
solve_poker_strategy()