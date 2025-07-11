def estimate_xpp_interpreter_runtime():
    """
    Calculates and prints the estimated worst-case runtime of an
    optimized C interpreter for the X++ language.
    """

    # --- Given Constants ---
    # Time costs in milliseconds (ms)
    T_READ_CHAR = 15
    T_RECOGNIZE_CHAR = 110

    # Worst-case scenario derived from the 366 character limit
    WORST_CASE_N = 121
    # The number '121' has 3 digits
    WORST_CASE_N_DIGITS = 3

    # --- Calculation ---

    # The fundamental cost is for processing a single character from the tape,
    # which involves reading it and then recognizing it.
    cost_per_char = T_READ_CHAR + T_RECOGNIZE_CHAR
    
    # 1. Time to parse the number of statements (n)
    # The first line of the input contains n. In the worst case, this is "121\n".
    # This line has 3 digits + 1 newline character = 4 characters total.
    chars_to_read_for_n = WORST_CASE_N_DIGITS + 1  # +1 for the newline
    time_read_n = chars_to_read_for_n * cost_per_char
    
    # 2. Time to execute the main loop for all statements
    # The loop runs n=121 times.
    # Each statement line (e.g., "X++\n") has 3 statement characters + 1 newline = 4 characters.
    chars_per_statement_line = 3 + 1
    cost_per_iteration = chars_per_statement_line * cost_per_char
    time_main_loop = WORST_CASE_N * cost_per_iteration
    
    # 3. Total Estimated Time
    # The total time is the sum of parsing n and executing the statements.
    # Printing the final integer is negligible (nanoseconds).
    total_time = time_read_n + time_main_loop
    
    # --- Output the detailed analysis ---
    
    print("--- Worst-Case Performance Analysis of Optimized X++ Interpreter ---")
    print("\n1. Time to parse the number of statements (n):")
    print(f"   - Worst-case n = {WORST_CASE_N}, which corresponds to the input line \"121\\n\".")
    print(f"   - This line has {chars_to_read_for_n} characters ({WORST_CASE_N_DIGITS} digits + 1 newline).")
    print(f"   - Cost per character = Read ({T_READ_CHAR} ms) + Recognize ({T_RECOGNIZE_CHAR} ms) = {cost_per_char} ms.")
    print(f"   - Time to read n = {chars_to_read_for_n} characters * {cost_per_char} ms/char = {time_read_n} ms.")

    print("\n2. Time to execute the main loop:")
    print(f"   - The loop runs {WORST_CASE_N} times.")
    print(f"   - Each statement line (e.g., \"X++\\n\") has {chars_per_statement_line} characters.")
    print(f"   - Cost per iteration = {chars_per_statement_line} characters * {cost_per_char} ms/char = {cost_per_iteration} ms.")
    print(f"   - Total loop time = {WORST_CASE_N} iterations * {cost_per_iteration} ms/iteration = {time_main_loop} ms.")

    print("\n3. Total Estimated Time:")
    print(f"   - Final Equation: {time_read_n} + {time_main_loop} = {total_time}")
    print(f"   - The estimated running time in the worst case is {total_time} ms.")

if __name__ == '__main__':
    estimate_xpp_interpreter_runtime()