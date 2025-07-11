import math

def calculate_and_print_time():
    """
    Calculates and prints the worst-case running time for the optimized X++ interpreter.
    """
    # Time costs in milliseconds (ms) as per the problem description
    T_READ_CHAR_MS = 15
    T_RECOGNIZE_CHAR_MS = 110
    T_COMPARE_CHAR_MS = 150

    # Total time for a single getchar() operation
    T_GETCHAR_TOTAL_MS = T_READ_CHAR_MS + T_RECOGNIZE_CHAR_MS

    # --- Step 1: Determine the worst-case (maximum) number of statements 'n' ---
    # The program on the tape is limited to 366 characters.
    # The format is: <number n><statement 1><statement 2>...<statement n>
    # Let 'd' be the number of digits in 'n'.
    # Each statement ("X++", "++X", etc.) is 3 characters long.
    # The total number of characters is: d + 3 * n <= 366
    # To find the maximum 'n', we test values.
    # For n = 121, d = 3. Total chars = 3 + 3 * 121 = 3 + 363 = 366. This fits.
    # For n = 122, d = 3. Total chars = 3 + 3 * 122 = 3 + 366 = 369. This is too large.
    # Therefore, the maximum number of statements is 121.
    n_max = 121
    digits_in_n_max = 3

    # --- Step 2: Calculate the time to read 'n' ---
    # Reading n=121 requires reading 3 digit characters.
    time_read_n_ms = digits_in_n_max * T_GETCHAR_TOTAL_MS

    # --- Step 3: Calculate the time to process one statement ---
    # The optimized interpreter reads 3 characters and performs 1 comparison per statement.
    # 1. Read first character.
    # 2. Read second character (the operator).
    # 3. Compare the second character to '+' or '-'.
    # 4. Read third character.
    time_per_statement_ms = (3 * T_GETCHAR_TOTAL_MS) + T_COMPARE_CHAR_MS

    # --- Step 4: Calculate the total time ---
    # Total time = (Time to read n) + (n_max * Time per statement)
    total_time_ms = time_read_n_ms + (n_max * time_per_statement_ms)

    # --- Step 5: Print the breakdown of the calculation ---
    print("Worst-case running time calculation:")
    print(f"Time for one getchar() operation = {T_READ_CHAR_MS} ms (read) + {T_RECOGNIZE_CHAR_MS} ms (recognize) = {T_GETCHAR_TOTAL_MS} ms")
    print(f"Maximum number of statements (n) fitting on a 366-char tape = {n_max}")
    print(f"Number of digits in n={n_max} is {digits_in_n_max}")
    print("\nFinal Equation:")
    print(f"T = (Time to read n) + (Number of statements) * (Time per statement)")
    print(f"T = ({digits_in_n_max} * {T_GETCHAR_TOTAL_MS}) + {n_max} * (3 * {T_GETCHAR_TOTAL_MS} + {T_COMPARE_CHAR_MS})")
    print(f"T = {time_read_n_ms} + {n_max} * ({3 * T_GETCHAR_TOTAL_MS} + {T_COMPARE_CHAR_MS})")
    print(f"T = {time_read_n_ms} + {n_max} * {time_per_statement_ms}")
    print(f"T = {time_read_n_ms} + {n_max * time_per_statement_ms}")
    print(f"T = {int(total_time_ms)} ms")

calculate_and_print_time()