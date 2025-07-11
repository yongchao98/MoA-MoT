import math

def solve_poker_strategy():
    """
    Calculates and explains the optimal poker strategy for the given scenario.
    """
    # Game parameters
    P = 10  # Pot size in dollars
    V = 6   # Hero's value combos (6 combos of AA)
    B = 1000 # Hero has $1000 behind, optimal bet is all-in

    # The optimal strategy is derived from making the villain indifferent to calling an all-in bet.
    # The GTO indifference equation is: N * (P + B) = V * B
    # where N is the number of bluff combos needed.
    # We solve for N to determine the required number of bluffs.
    
    print("To find the optimal strategy, we first determine the optimal bet size for our value hand (AA).")
    print("To maximize value, we should bet as large as possible, which is all-in ($1000).")
    print("Next, we must balance this bet with bluffs to get called. We use the GTO indifference formula to find the number of bluff combos (N) needed.")
    print("\nIndifference Equation: N * (Pot + Bet) = ValueCombos * Bet")
    
    # Printing each number in the final equation
    print(f"N * ({P} + {B}) = {V} * {B}")
    
    term1 = P + B
    term2 = V * B
    print(f"N * {term1} = {term2}")
    
    N = term2 / term1
    print(f"N = {term2} / {term1}")
    print(f"N = {N:.2f}")

    print("\nThis means Hero needs to bluff with approximately 5.94 combos.")
    print("Since Hero has 60 available bluff combos (QQ-33), this is a sustainable strategy.")
    print("Therefore, the optimal strategy for AA is to bet all-in to maximize value.")
    
    # Final strategy for AA, rounded to nearest even numbers per instructions
    action = "BET"
    sizing = B
    frequency_percent = 100
    
    final_sizing = int(round(sizing / 2.0) * 2)
    final_frequency = int(round(frequency_percent / 2.0) * 2)
    
    final_answer_string = f"{action} {final_sizing} {final_frequency}%"
    
    print("\nYour optimal strategy with AA is:")
    print(final_answer_string)
    
    # Output the final answer in the requested format
    print(f"\n<<<{final_answer_string}>>>")

solve_poker_strategy()