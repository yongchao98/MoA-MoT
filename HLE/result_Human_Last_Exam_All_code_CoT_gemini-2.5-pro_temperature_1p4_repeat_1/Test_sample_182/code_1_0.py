import math

def calculate_and_print_worst_case_time():
    """
    Calculates and prints the step-by-step estimation for the
    worst-case running time of the optimized C interpreter.
    """
    # Step 1: Define performance costs in milliseconds (ms)
    time_read_char = 15
    time_recognize_char = 110
    time_compare_char = 150

    # Combined time to read and recognize a single character
    time_get_char = time_read_char + time_recognize_char

    # Step 2: Determine the worst-case number of statements (n)
    # The total number of characters on the tape is limited to 366.
    # Program length = (number of digits in n) + 3 * n
    # By solving d + 3n <= 366, we find the maximum n is 121.
    worst_case_n = 121
    num_digits_in_n = len(str(worst_case_n))

    # Step 3: Calculate the cost for each part of the execution

    # Time to read the number 'n' (e.g., "121")
    time_to_read_n = num_digits_in_n * time_get_char

    # Time to process a single statement with the optimized algorithm:
    # Read 1st char + Read 2nd char + Compare 2nd char + Read 3rd char
    time_per_statement = time_get_char + time_get_char + time_compare_char + time_get_char

    # Time to process all statements
    total_time_for_statements = worst_case_n * time_per_statement

    # Total worst-case execution time
    total_worst_case_time = time_to_read_n + total_time_for_statements

    # Step 4: Print the breakdown of the final calculation
    print("Worst-Case Running Time Estimation:")
    print("---------------------------------")
    print(f"Maximum number of statements (n): {worst_case_n}")
    print(f"Number of digits in n (d): {num_digits_in_n}")
    print(f"Time to get a single character: {time_get_char} ms")
    print(f"Time to compare two characters: {time_compare_char} ms")
    print("\nFinal Equation:")
    print("Total Time = (d * Time_get_char) + (n * (3 * Time_get_char + Time_compare_char))")
    print(f"Total Time = ({num_digits_in_n} * {time_get_char}) + ({worst_case_n} * (3 * {time_get_char} + {time_compare_char}))")
    # Note: I am expanding the (3*T+T_comp) term below for clarity, as per the optimized logic.
    print(f"Total Time = ({num_digits_in_n} * {time_get_char}) + ({worst_case_n} * ({time_get_char} + {time_get_char} + {time_compare_char} + {time_get_char}))")
    print(f"Total Time = ({time_to_read_n}) + ({worst_case_n} * {time_per_statement})")
    print(f"Total Time = {time_to_read_n} + {total_time_for_statements}")
    print(f"Total Time = {total_worst_case_time} ms")

calculate_and_print_worst_case_time()