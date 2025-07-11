import math

# Plan:
# 1. Define the time costs for elemental operations in milliseconds.
# 2. Determine the worst-case number of statements (n) that fits the 366 character limit.
# 3. Based on an optimized C interpreter logic, count the number of expensive operations (getchar, char comparison).
# 4. Calculate the total time and print the breakdown.

# --- Step 1: Define time costs ---
T_read_char_ms = 15
T_recognize_char_ms = 110
T_compare_char_ms = 150
T_int_op_ns = 10 # Negligible
T_print_char_ns = 20 # Negligible

# The total time for a single getchar() operation is reading + recognition.
T_getchar_ms = T_read_char_ms + T_recognize_char_ms

# --- Step 2: Determine worst-case n ---
# Program size = (digits_in_n + 1) + (n * 4) <= 366
# To maximize n, we test for a small number of digits in n.
# For n=90 (2 digits): (2 + 1) + (90 * 4) = 3 + 360 = 363. This fits.
n_max = 90
digits_in_n_max = 2

# --- Step 3: Count operations for optimized C code ---
# Optimized Logic for n=90:
# a. Read n: Read '9', '0', and newline. This is 3 getchar() calls.
# b. Loop 90 times:
#    - Read char 1 (discard)
#    - Read char 2 (use for logic)
#    - Compare char 2 to '+'
#    - Read char 3 (discard)
#    - Read newline (discard)
#    This is 4 getchar() calls and 1 char comparison per loop.

# Total getchar() calls:
getchar_for_n = digits_in_n_max + 1  # For '9', '0', and newline
getchar_for_statements = n_max * 4   # 4 reads per statement line
total_getchar_calls = getchar_for_n + getchar_for_statements

# Total character comparisons:
# One comparison per statement to check if the middle char is '+' or '-'.
total_char_comparisons = n_max

# --- Step 4: Calculate and Print Total Time ---
time_from_getchar = total_getchar_calls * T_getchar_ms
time_from_comparisons = total_char_comparisons * T_compare_char_ms
total_time_ms = time_from_getchar + time_from_comparisons

print("Worst-Case Running Time Estimation for Optimized C Interpreter")
print("-" * 60)
print(f"Maximum number of statements (n) fitting in 366 characters: {n_max}")
print(f"Cost per getchar() = {T_read_char_ms}ms (read) + {T_recognize_char_ms}ms (recognize) = {T_getchar_ms} ms")
print(f"Cost per character comparison = {T_compare_char_ms} ms")
print("-" * 60)
print("Calculation Breakdown:")
print(f"Total getchar() calls = (for n='{n_max}') + (for {n_max} statements)")
print(f"                      = ({digits_in_n_max} + 1) + ({n_max} * 4) = {total_getchar_calls} calls")
print(f"Total character comparisons = {n_max} statements * 1 comparison/statement = {total_char_comparisons} comparisons")
print("")
print("Final Time Equation:")
print(f"Total Time = (Total getchar() calls * Time per getchar()) + (Total comparisons * Time per comparison)")
print(f"Total Time = ({total_getchar_calls} * {T_getchar_ms}) + ({total_char_comparisons} * {T_compare_char_ms})")
print(f"Total Time = ({time_from_getchar}) + ({time_from_comparisons})")
print(f"Total Time = {total_time_ms} ms")
