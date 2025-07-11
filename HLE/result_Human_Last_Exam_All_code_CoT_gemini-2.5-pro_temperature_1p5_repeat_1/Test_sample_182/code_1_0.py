# Step 1: Define the performance costs in milliseconds from the problem description.
time_read_from_tape_ms = 15
time_recognize_symbol_ms = 110
time_compare_two_chars_ms = 150

# The cost of the getchar() operation includes both reading and recognition.
time_per_getchar_ms = time_read_from_tape_ms + time_recognize_symbol_ms

# Step 2: Determine the worst-case scenario (maximum number of statements, n).
# The total number of characters on the tape cannot exceed 366.
# The size of the program is: number_of_digits(n) + 1 (for newline) + n * (3 + 1).
# Let d be the number of digits in n. The inequality is: d + 1 + 4*n <= 366.
# By testing values, we find the maximum n occurs when n has 2 digits (d=2).
# 2 + 1 + 4*n <= 366  =>  4*n <= 363  =>  n <= 90.75.
# Thus, the maximum integer n is 90.
worst_case_n = 90
digits_in_n = 2

# Step 3: Calculate the total execution time for the optimized interpreter.

# Part A: Time to read the number of statements, n.
# For n=90, the interpreter must read 2 digits and 1 newline character.
chars_to_read_n = digits_in_n + 1
time_to_read_n_ms = chars_to_read_n * time_per_getchar_ms

# Part B: Time to process all n statements in the main loop.
# The optimized algorithm for each statement is:
# 1. Read 3 statement characters + 1 newline character = 4 getchar() calls.
# 2. Perform 1 character comparison (on the second character, which is always '+' or '-')
time_per_statement_ms = (4 * time_per_getchar_ms) + time_compare_two_chars_ms
total_time_for_loop_ms = worst_case_n * time_per_statement_ms

# The total estimated time is the sum of Part A and Part B.
# Printing the final integer result is considered negligible (20 ns per character).
# The final equation is: total_time = time_to_read_n_ms + total_time_for_loop_ms
# The numbers in the equation are time_to_read_n_ms and total_time_for_loop_ms.
total_estimated_time_ms = time_to_read_n_ms + total_time_for_loop_ms

# Print the final result, which is the integer value of T.
print(int(total_estimated_time_ms))