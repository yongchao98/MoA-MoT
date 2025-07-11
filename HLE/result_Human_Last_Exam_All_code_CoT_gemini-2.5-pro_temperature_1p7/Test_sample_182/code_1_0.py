import math

def calculate_worst_case_runtime():
    """
    Calculates the estimated worst-case running time for the optimized C interpreter.
    """
    # Timing constants provided in the problem, in milliseconds (ms)
    T_READ_CHAR = 15
    T_RECOGNIZE_DIGIT = 110
    T_COMPARE_CHAR = 150

    # Step 1: Determine the worst-case number of statements (n)
    # The tape limit is 366 characters.
    # Total characters = D (digits of n) + 1 (newline) + n * 4 (statement + newline)
    # Testing for n=90 (D=2): 2 + 1 + 90 * 4 = 363 <= 366. This is valid.
    # Testing for n=91 (D=2): 2 + 1 + 91 * 4 = 367 > 366. This is invalid.
    # Thus, the worst-case n is 90.
    N_WORST = 90
    D_WORST = 2  # Number of digits in 90

    # Step 2: Calculate the time to read the number 'n' from the tape
    # The logic is to read character by character until a newline is found.
    # For each digit read, we perform a read, a comparison (to check for newline),
    # and a recognition (to convert character to digit).
    time_per_digit = T_READ_CHAR + T_COMPARE_CHAR + T_RECOGNIZE_DIGIT
    
    # For the final newline character, we just do a read and a comparison.
    time_for_newline_read = T_READ_CHAR + T_COMPARE_CHAR
    
    time_read_n = (D_WORST * time_per_digit) + time_for_newline_read

    # Step 3: Calculate the time for the main loop processing n statements
    # The optimized plan for each statement is:
    # 1. getchar() -> Read 1st char
    # 2. c = getchar() -> Read 2nd char
    # 3. if (c == '+') -> Compare 2nd char
    # 4. getchar() -> Read 3rd char
    # 5. getchar() -> Read newline
    # This amounts to 4 character reads and 1 character comparison per statement.
    time_per_statement = (4 * T_READ_CHAR) + (1 * T_COMPARE_CHAR)
    time_main_loop = N_WORST * time_per_statement

    # Step 4: Calculate the total estimated time
    total_time = time_read_n + time_main_loop
    
    # --- Output the step-by-step calculation ---

    print("--- Worst-Case Runtime Estimation for Optimized C Interpreter ---")

    print("\n[A] Time to read n (for n=90):")
    print(f"Time per digit = Read_Char + Compare_Char + Recognize_Digit")
    print(f"               = {T_READ_CHAR} + {T_COMPARE_CHAR} + {T_RECOGNIZE_DIGIT} = {time_per_digit} ms")
    print(f"Time for final newline = Read_Char + Compare_Char")
    print(f"                       = {T_READ_CHAR} + {T_COMPARE_CHAR} = {time_for_newline_read} ms")
    print(f"Total time to read n = (Num_Digits * Time_per_digit) + Time_for_final_newline")
    print(f"                     = ({D_WORST} * {time_per_digit}) + {time_for_newline_read} = {time_read_n} ms")

    print("\n[B] Time for main loop (for 90 statements):")
    print(f"Time per statement = (4 * Read_Char) + (1 * Compare_Char)")
    print(f"                   = (4 * {T_READ_CHAR}) + (1 * {T_COMPARE_CHAR}) = {time_per_statement} ms")
    print(f"Total time for loop = Num_Statements * Time_per_statement")
    print(f"                    = {N_WORST} * {time_per_statement} = {time_main_loop} ms")

    print("\n[C] Total Estimated Time:")
    print(f"Total Time = Time_to_read_n + Time_for_main_loop")
    print(f"           = {time_read_n} + {time_main_loop} = {int(total_time)} ms")

calculate_worst_case_runtime()