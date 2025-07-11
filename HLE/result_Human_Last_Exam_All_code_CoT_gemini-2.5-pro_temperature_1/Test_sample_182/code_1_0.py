def calculate_worst_case_runtime():
    """
    Calculates the worst-case running time for the optimized X++ interpreter.
    """
    # --- Given Performance Costs (in milliseconds) ---
    T_read_char = 15
    T_recognize_char = 110
    T_compare_char = 150

    # --- Derived Costs ---
    # Cost to read a character and recognize what it is
    T_read_and_recognize = T_read_char + T_recognize_char
    # Cost to read a character that we are just discarding (e.g., end of line)
    T_read_and_discard = T_read_char

    # --- Worst-Case Scenario Analysis ---
    # 1. Determine the maximum number of statements (n)
    # Total chars on tape <= 366.
    # Program chars = (digits_in_n + 1 for newline) + (n * 4 for each statement "X++\n")
    # For n=90 (2 digits), chars = (2 + 1) + (90 * 4) = 3 + 360 = 363. This fits.
    # For n=91 (2 digits), chars = (2 + 1) + (91 * 4) = 3 + 364 = 367. This is too large.
    # So, the maximum number of statements is 90.
    n_max = 90
    digits_in_n_max = 2

    # 2. Calculate time to read n
    # The C code reads n by reading and recognizing each digit.
    # For n=90, it reads '9' and '0'.
    # It then consumes the trailing newline character.
    T_read_n_value = digits_in_n_max * T_read_and_recognize
    T_consume_newline_after_n = T_read_and_discard
    T_setup = T_read_n_value + T_consume_newline_after_n
    
    print("--- Worst-Case Runtime Calculation ---")
    print("\nStep 1: Calculate time to read the number of statements (n=90)")
    print(f"Time to read and recognize {digits_in_n_max} digits: {digits_in_n_max} * ({T_read_char} + {T_recognize_char}) = {T_read_n_value} ms")
    print(f"Time to consume newline after n: {T_consume_newline_after_n} ms")
    print(f"Total setup time: {T_read_n_value} + {T_consume_newline_after_n} = {T_setup} ms")

    # 3. Calculate time for one worst-case statement
    # The optimized logic reads 2 chars (c1, c2), then checks 'if (c1=='+' || c2=='+')'.
    # A worst-case statement like "X--" makes both comparisons evaluate to false.
    T_read_2_chars = 2 * T_read_and_recognize
    T_compare_2_chars = 2 * T_compare_char
    # Discard the rest of the line (1 remaining char + 1 newline char)
    T_discard_rest_of_line = 2 * T_read_and_discard
    T_statement_worst = T_read_2_chars + T_compare_2_chars + T_discard_rest_of_line

    print("\nStep 2: Calculate time for a single worst-case statement (e.g., 'X--')")
    print(f"Time to read and recognize first 2 chars: 2 * ({T_read_char} + {T_recognize_char}) = {T_read_2_chars} ms")
    print(f"Time for 2 character comparisons in the 'if': 2 * {T_compare_char} = {T_compare_2_chars} ms")
    print(f"Time to discard remaining 2 chars on line: 2 * {T_read_and_discard} = {T_discard_rest_of_line} ms")
    print(f"Total time per worst-case statement: {T_read_2_chars} + {T_compare_2_chars} + {T_discard_rest_of_line} = {T_statement_worst} ms")

    # 4. Calculate total worst-case time
    T_loop = n_max * T_statement_worst
    T_total = T_setup + T_loop

    print("\nStep 3: Calculate total worst-case time for the entire program")
    print(f"Total time for the loop ({n_max} statements): {n_max} * {T_statement_worst} = {T_loop} ms")
    print(f"\nFinal Equation for Total Time:")
    print(f"Total Time = (Time to read n) + (n * Time per statement)")
    print(f"Total Time = {T_setup} + {n_max} * {T_statement_worst} = {T_total} ms")
    
    print(f"\nEstimated running time in the worst case (rounded):")
    print(f"{round(T_total)} ms")

calculate_worst_case_runtime()