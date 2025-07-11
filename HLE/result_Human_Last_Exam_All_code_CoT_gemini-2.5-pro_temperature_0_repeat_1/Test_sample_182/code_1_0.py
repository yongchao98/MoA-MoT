def estimate_worst_case_runtime():
    """
    Calculates the estimated worst-case running time for the optimized X++ interpreter.
    """
    # --- Given Costs (in milliseconds) ---
    read_char_time = 15
    recognize_char_time = 110
    # Integer and printing operations are in nanoseconds, so they are negligible.

    # --- Step 1: Determine the maximum number of statements (n) ---
    # The total program length is at most 366 characters.
    # A program consists of a line for n (d digits + 1 newline) and n statement lines (3 chars + 1 newline).
    # Total characters = d + 1 + n * 4.
    # We need to find the max n such that d + 1 + 4*n <= 366.
    # Let's test for a 2-digit n (d=2): 2 + 1 + 4*n <= 366 => 4n <= 363 => n <= 90.75.
    # So, the maximum possible value for n is 90.
    max_n = 90
    digits_in_n = 2 # "90" has 2 digits

    # --- Step 2: Calculate the time to read n ---
    # The worst case for n=90 means the first line is "90" followed by a newline.
    # We must read and recognize '9' and '0'. We only need to read the newline.
    time_to_read_n = (read_char_time + recognize_char_time) * digits_in_n + read_char_time * 1

    # --- Step 3: Calculate the time for one worst-case statement ---
    # A worst-case statement is 'X++' or 'X--'. Let's analyze 'X++'.
    # 1. Read 'X': 1 read (15ms) + 1 recognition (110ms) = 125ms
    # 2. Read '+': 1 read (15ms) + 1 recognition (110ms) = 125ms
    #    (At this point, the operation is identified)
    # 3. Consume the rest of the line (the second '+' and the newline)
    #    Read '+': 1 read (15ms)
    #    Read newline: 1 read (15ms)
    # Total time for one statement:
    time_per_statement = (read_char_time + recognize_char_time) * 2 + read_char_time * 2

    # --- Step 4: Calculate the total time ---
    total_time = time_to_read_n + max_n * time_per_statement

    # --- Print the final equation and result ---
    print("Worst-case running time calculation:")
    print(f"Time to read n (for n={max_n}): {time_to_read_n} ms")
    print(f"Number of statements (n): {max_n}")
    print(f"Time per worst-case statement: {time_per_statement} ms")
    print("\nFinal Equation:")
    print(f"Total Time = (Time to read n) + n * (Time per statement)")
    print(f"Total Time = {time_to_read_n} + {max_n} * {time_per_statement}")
    final_result = int(time_to_read_n + max_n * time_per_statement)
    print(f"Total Time = {final_result} ms")

estimate_worst_case_runtime()