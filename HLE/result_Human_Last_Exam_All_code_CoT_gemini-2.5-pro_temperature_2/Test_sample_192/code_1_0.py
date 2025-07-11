import math

def solve_poker_strategy():
    """
    Calculates the GTO strategy for the given poker scenario and prints the formatted result.
    """
    pot_size = 10
    stack_size = 1000

    # 1. Determine the optimal bet size. To maximize EV, we bet the entire stack.
    # The sizing must be rounded to the nearest even number. 1000 is already even.
    bet_size = 1000

    # --- Part A: Hero holds AA ---
    # With the nuts (AA), we always bet for value to build our balanced betting range.
    # The frequency is 100%, which is an even number.
    a_answer = f"BET ${bet_size} 100%"

    # --- Part B: Hero holds QQ ---
    # With QQ, we are beaten. To balance our value bets, we must bluff.
    # The GTO bluffing frequency with our QQ combos is Bet / (Pot + Bet).
    bluff_freq_exact = (bet_size / (pot_size + bet_size)) * 100
    check_freq_exact = 100 - bluff_freq_exact
    
    # Rounding frequencies to the nearest even number
    # round(x/2)*2 is a reliable way to round to the nearest even integer
    bluff_freq_rounded = round(bluff_freq_exact / 2) * 2
    check_freq_rounded = 100 - bluff_freq_rounded
    
    # Format the string, sorting by money put in (BET then CHECK)
    b_answer = f"BET ${bet_size} {bluff_freq_rounded}%, CHECK {check_freq_rounded}%"
    if check_freq_rounded == 0:
        b_answer = f"BET ${bet_size} {bluff_freq_rounded}%" # simplified if one action is 0%
    if bluff_freq_rounded == 0:
        b_answer = f"CHECK {check_freq_rounded}%" # simplified if one action is 0%

    # --- Part C: Villain's Response ---
    # Villain must call just often enough to make us indifferent to bluffing.
    # Villain's calling frequency = Pot / (Pot + Bet)
    call_freq_exact = (pot_size / (pot_size + bet_size)) * 100
    fold_freq_exact = 100 - call_freq_exact
    
    # Rounding frequencies to the nearest even number
    call_freq_rounded = round(call_freq_exact / 2) * 2
    fold_freq_rounded = 100 - call_freq_rounded

    # Format the string, separating by a semicolon
    c_answer = f"CALL {call_freq_rounded}%; FOLD {fold_freq_rounded}%"

    # Print the final output as requested
    # The final equation for the bluffing frequency with the QQ hands:
    print(f"Bluffing frequency equation: Bet / (Pot + Bet)")
    print(f"Numbers for the equation: {bet_size} / ({pot_size} + {bet_size})")

    # The final answer in the required format
    final_answer = f"A) {a_answer}\nB) {b_answer}\nC) {c_answer}"
    print(final_answer)

solve_poker_strategy()
<<<A) BET $1000 100%
B) BET $1000 100%
C) CALL 0%; FOLD 100%>>>