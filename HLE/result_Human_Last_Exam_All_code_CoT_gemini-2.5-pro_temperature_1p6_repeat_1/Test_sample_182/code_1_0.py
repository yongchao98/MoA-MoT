import math

def calculate_worst_case_runtime():
    """
    Calculates and prints the step-by-step estimation of the worst-case
    running time for the optimized X++ interpreter.
    """
    # --- Step 1: Define cost constants from the problem description ---
    t_read_char_ms = 15
    t_recognize_ms = 110
    
    # --- Step 2: Determine worst-case parameters ---
    # The maximum number of statements (n) that fits on a 366-character tape.
    # Total chars = (digits in n) + 1 (newline) + n * (3 chars/statement + 1 newline)
    # By testing, n=90 is the maximum. (2 + 1 + 90 * 4 = 363 <= 366)
    n_max = 90
    digits_in_n_max = 2

    print("### Worst-Case Running Time Estimation ###")
    print("\n# This script calculates the total time based on the cost model.")
    print("\n# Cost per operation (in ms):")
    print(f"T_read_char = {t_read_char_ms}")
    print(f"T_recognize = {t_recognize_ms}")

    # --- Step 3: Calculate the time to read the number of statements (n) ---
    # For n=90, we read '9' and '0'. This requires 2 reads and 2 recognitions.
    # We also read one end-of-line character.
    time_read_n = digits_in_n_max * (t_read_char_ms + t_recognize_ms) + t_read_char_ms

    print("\n# 1. Time to read the number of statements (n_max = 90):")
    print(f"   Equation: n_digits * (T_read_char + T_recognize) + T_read_char_for_newline")
    print(f"   Calculation: {digits_in_n_max} * ({t_read_char_ms} + {t_recognize_ms}) + {t_read_char_ms}")
    print(f"   Result: {time_read_n} ms")

    # --- Step 4: Calculate the time for one worst-case statement ---
    # The worst-case statement is "X++" or "X--".
    # a) Read 'X' (15ms) + Recognize 'X' (110ms)
    # b) Read '+' (15ms) + Recognize '+' (110ms)
    # c) The operation is now known. To advance the tape, read the remaining
    #    two characters on the line (the third character and the newline).
    #    This takes 2 more character reads.
    time_per_statement_worst = 2 * (t_read_char_ms + t_recognize_ms) + 2 * t_read_char_ms

    print("\n# 2. Time to process one worst-case statement (e.g., 'X++'):")
    print("   Equation: 2*(T_read_char + T_recognize) + 2*T_read_char")
    print(f"   Calculation: 2*({t_read_char_ms} + {t_recognize_ms}) + 2*{t_read_char_ms}")
    print(f"   Result: {time_per_statement_worst} ms")

    # --- Step 5: Calculate total execution time ---
    # Total time = (time to read n) + n_max * (time for one worst-case statement)
    time_loop = n_max * time_per_statement_worst
    total_time = time_read_n + time_loop

    print("\n# 3. Total time for executing the loop of n_max statements:")
    print("   Equation: n_max * Time_per_statement")
    print(f"   Calculation: {n_max} * {time_per_statement_worst}")
    print(f"   Result: {time_loop} ms")

    print("\n# 4. Total worst-case running time:")
    print("   Equation: Time_read_n + Time_loop")
    print(f"   Calculation: {time_read_n} + {time_loop}")
    print(f"   Result: {total_time} ms")
    
    # Final answer as per instructions
    final_answer = round(total_time)
    print(f"\nFinal rounded answer is {final_answer} ms.")

calculate_worst_case_runtime()