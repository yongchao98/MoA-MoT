def calculate_worst_case_time():
    """
    Calculates the worst-case running time for the optimized X++ interpreter.
    """
    # Cost constants in milliseconds (ms)
    T_read_char = 15
    T_recog_char = 110

    # Step 1: Determine the maximum number of statements (n)
    # A program's character count is len(str(n)) + 1 (newline) + n * (3-char stmt + 1 newline).
    # We test values for n to find the maximum that fits under 366 characters.
    # For n=90, len("90") + 1 + 90 * 4 = 2 + 1 + 360 = 363 <= 366 (This works)
    # For n=91, len("91") + 1 + 91 * 4 = 2 + 1 + 364 = 367 > 366 (This fails)
    # Therefore, the maximum number of statements in the worst case is 90.
    n_max = 90
    n_digits = 2  # The number "90" has 2 digits.

    # Step 2: Calculate the time to read the number 'n' from the tape.
    # This involves reading and recognizing each digit, plus reading the final newline.
    # Time = (T_read + T_recog for '9') + (T_read + T_recog for '0') + (T_read for '\n')
    time_read_n = n_digits * (T_read_char + T_recog_char) + T_read_char

    # Step 3: Calculate the time to process the single worst-case statement.
    # The worst-case statement is one like 'X++' or 'X--', where the operation
    # symbol is not the first character.
    # Read 'X', Recognize 'X' -> The interpreter must check if it's '+' or '-'.
    # Read '+', Recognize '+' -> The operation is now identified.
    # Read '+' (2nd one)      -> No recognition needed, op is already known.
    # Read '\n'               -> No recognition needed.
    time_per_worst_statement = (T_read_char + T_recog_char) + \
                               (T_read_char + T_recog_char) + \
                               T_read_char + \
                               T_read_char

    # Step 4: Calculate the total time for processing all statements.
    time_all_statements = n_max * time_per_worst_statement

    # Step 5: Calculate the total worst-case time.
    total_time = time_read_n + time_all_statements

    # Step 6: Print the detailed breakdown of the calculation.
    print(f"Worst-case analysis for the optimized X++ interpreter:")
    print(f"Maximum number of statements (n) fitting on the tape: {n_max}")
    print("-" * 50)
    print("1. Time to read n = '90':")
    print(f"   Reading and recognizing {n_digits} digits = {n_digits} * ({T_read_char}ms + {T_recog_char}ms) = {n_digits * (T_read_char + T_recog_char)} ms")
    print(f"   Reading the newline character = {T_read_char} ms")
    print(f"   Subtotal for reading n: {time_read_n} ms")
    print("-" * 50)
    print(f"2. Time to process {n_max} worst-case statements (e.g., 'X++'):")
    print(f"   Time per statement = ({T_read_char}+{T_recog_char}) + ({T_read_char}+{T_recog_char}) + {T_read_char} + {T_read_char} = {time_per_worst_statement} ms")
    print(f"   Subtotal for {n_max} statements = {n_max} * {time_per_worst_statement} ms = {time_all_statements} ms")
    print("-" * 50)
    print("Final Calculation:")
    print(f"Total Time = (Time to read n) + (Time for all statements)")
    print(f"Total Time = {time_read_n} + {time_all_statements} = {int(total_time)}")
    print(f"\nThe estimated running time in the worst case is {int(total_time)} ms.")

calculate_worst_case_time()