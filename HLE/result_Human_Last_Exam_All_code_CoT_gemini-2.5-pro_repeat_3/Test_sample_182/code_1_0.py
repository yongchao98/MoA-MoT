import math

def calculate_worst_case_time():
    """
    Calculates the estimated worst-case running time for the optimized X++ interpreter.
    """
    # Performance constants from the problem description (in milliseconds)
    T_READ_CHAR = 15
    T_RECOG_CHAR = 110
    T_COMPARE_CHAR = 150

    # Combined time to read and recognize a single character
    T_GETCHAR = T_READ_CHAR + T_RECOG_CHAR

    # --- Step 1: Determine the worst-case number of statements (n) ---
    # The total number of characters on the tape is at most 366.
    # Program format: number `n` followed by `n` statements of 3 characters each.
    # Total characters = len(str(n)) + 3 * n.
    # We need to find the max n such that len(str(n)) + 3 * n <= 366.
    # By inspection:
    # for n = 121, len("121") + 3 * 121 = 3 + 363 = 366. This is the maximum.
    max_n = 121
    n_str = str(max_n)
    len_n_str = len(n_str)

    # --- Step 2: Calculate the time to read n (T_read_n) ---
    # To read n="121", the program reads 3 digits.
    # For each digit, it performs one getchar() and two comparisons (c >= '0' and c <= '9').
    cost_per_digit = T_GETCHAR + 2 * T_COMPARE_CHAR
    time_for_digits = len_n_str * cost_per_digit
    
    # An additional getchar() is needed to read the end-of-line after the number.
    time_for_eol_after_n = T_GETCHAR
    
    T_read_n = time_for_digits + time_for_eol_after_n

    # --- Step 3: Calculate the time for the worst-case statement ---
    # The optimized interpreter reads 3 characters per statement. It stops checking for
    # an operator once it has found one.
    # Worst case statement: "X++" or "X--", where the operator appears at the 2nd position.
    
    # Cost for Char 1 ('X'): 1 getchar() + 1 compare ('+') + 1 compare ('-')
    cost_char1 = T_GETCHAR + T_COMPARE_CHAR + T_COMPARE_CHAR
    
    # Cost for Char 2 ('+'): 1 getchar() + 1 compare ('+') -> operator found
    cost_char2 = T_GETCHAR + T_COMPARE_CHAR
    
    # Cost for Char 3 ('+'): 1 getchar() -> operator already found, no comparisons needed
    cost_char3 = T_GETCHAR
    
    T_statement_worst = cost_char1 + cost_char2 + cost_char3

    # --- Step 4: Calculate total program execution time ---
    # The main loop runs 'max_n' times.
    T_loop = max_n * T_statement_worst
    
    # Total time is the sum of reading n and executing the loop. Printing time is negligible.
    T_total = T_read_n + T_loop
    
    # The result should be an integer.
    final_answer = int(round(T_total))

    # --- Step 5: Print the results as requested ---
    print("Worst-case analysis of the optimized interpreter:")
    print(f"Max statements (n): {max_n}")
    print(f"Time to read n: {T_read_n} ms")
    print(f"Time for one worst-case statement: {T_statement_worst} ms")
    print(f"Total time for statement loop: {T_loop} ms")
    print("-" * 30)
    print("Final equation for total time (T):")
    print(f"T = (Time_to_read_n) + (n * Time_per_worst_statement)")
    print(f"T = ({len_n_str} * ({T_GETCHAR} + 2 * {T_COMPARE_CHAR}) + {T_GETCHAR}) + {max_n} * (({T_GETCHAR} + 2 * {T_COMPARE_CHAR}) + ({T_GETCHAR} + {T_COMPARE_CHAR}) + {T_GETCHAR})")
    print(f"T = {T_read_n} + {max_n} * {T_statement_worst}")
    print(f"T = {T_read_n} + {T_loop}")
    print(f"T = {final_answer} ms")


if __name__ == '__main__':
    calculate_worst_case_time()