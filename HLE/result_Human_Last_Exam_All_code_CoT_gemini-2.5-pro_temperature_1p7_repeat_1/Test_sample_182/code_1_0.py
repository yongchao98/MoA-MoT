import math

def calculate_worst_case_runtime():
    """
    Calculates the estimated worst-case runtime for the optimized X++ interpreter.
    """
    # --- Step 1: Define the time constants from the problem description (in milliseconds) ---
    time_read_char_ms = 15
    time_recognize_digit_ms = 110
    time_compare_char_ms = 150
    # Integer and print operations are in nanoseconds, which are negligible compared to millisecond operations.

    # --- Step 2: Find the maximum number of statements (n) that fit on the tape ---
    max_tape_chars = 366
    max_n = 0
    # Iterate to find the largest n where num_digits(n) + 3*n <= 366
    for n_candidate in range(1, 200): # A sufficiently large range to check
        num_digits = len(str(n_candidate))
        total_chars_needed = num_digits + (3 * n_candidate)
        if total_chars_needed <= max_tape_chars:
            max_n = n_candidate
        else:
            break # Stop when we exceed the tape limit

    num_digits_in_max_n = len(str(max_n))

    # --- Step 3: Calculate the time to read n in the worst case (i.e., for max_n) ---
    time_to_read_n = num_digits_in_max_n * (time_read_char_ms + time_recognize_digit_ms)

    # --- Step 4: Calculate the time to process one statement in the worst case ---
    # The worst case is a statement like "X--" which forces the maximum number of comparisons.
    # 1. Read 'X'                           (15 ms)
    # 2. Compare 'X' with '+' (fails)        (150 ms)
    # 3. Compare 'X' with '-' (fails)        (150 ms)
    # 4. Read the second char, '-'           (15 ms)
    # 5. Compare '-' with '+' (fails)        (150 ms)
    # 6. Operation is identified as decrement.
    # 7. Read the third char, '-', to clear the line (15 ms)
    time_per_worst_statement = (
        time_read_char_ms +     # Read 1st char
        time_compare_char_ms +  # Compare with '+'
        time_compare_char_ms +  # Compare with '-'
        time_read_char_ms +     # Read 2nd char
        time_compare_char_ms +  # Compare with '+'
        time_read_char_ms       # Read 3rd char
    )

    # --- Step 5: Calculate the total worst-case time ---
    total_time = time_to_read_n + (max_n * time_per_worst_statement)

    # --- Step 6: Print the detailed calculation and the final result ---
    print("Worst-Case Runtime Estimation:")
    print("-" * 30)

    print(f"1. Maximum number of statements (n) fitting on a {max_tape_chars} character tape:")
    print(f"   The largest 'n' where len(str(n)) + 3*n <= {max_tape_chars} is: {max_n}")
    print()

    print(f"2. Time to read the number n = {max_n} ({num_digits_in_max_n} digits):")
    print(f"   Equation: {num_digits_in_max_n} * (read_char_time + recognize_digit_time)")
    print(f"   Calculation: {num_digits_in_max_n} * ({time_read_char_ms} ms + {time_recognize_digit_ms} ms) = {time_to_read_n} ms")
    print()

    print("3. Time for the single worst-case statement (e.g., 'X--'):")
    print(f"   Equation: read + cmp + cmp + read + cmp + read")
    print(f"   Calculation: {time_read_char_ms} + {time_compare_char_ms} + {time_compare_char_ms} + {time_read_char_ms} + {time_compare_char_ms} + {time_read_char_ms} = {time_per_worst_statement} ms")
    print()

    print("4. Total worst-case time:")
    print(f"   Equation: time_to_read_n + (n * time_per_statement)")
    print(f"   Calculation: {time_to_read_n} + ({max_n} * {time_per_worst_statement}) = {total_time} ms")
    print("-" * 30)
    
    # Final result, rounded to the nearest integer as requested
    final_answer = round(total_time)
    print(f"Final Answer: {final_answer} ms")


if __name__ == "__main__":
    calculate_worst_case_runtime()