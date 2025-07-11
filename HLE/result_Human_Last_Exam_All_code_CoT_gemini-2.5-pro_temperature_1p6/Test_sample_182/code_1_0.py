import math

# Step 1: Define the performance costs in milliseconds from the problem description.
# Reading from tape takes 15ms, and recognizing a symbol takes 110ms.
T_READ_CHAR = 15 + 110
# Comparing two characters takes 150ms.
T_COMPARE_CHAR = 150

# Step 2: Define the parameters for the worst-case scenario.
# As explained in the plan, the maximum number of statements 'n' that fits
# the 366-character tape limit is 90.
N_MAX = 90
# The number of digits in n=90 is 2.
NUM_DIGITS_N = 2

# Step 3: Calculate the total time for the optimized interpreter.

# Time to read the number 'n' (e.g., "90\n"). This involves reading
# each digit plus the newline character.
time_read_n = (NUM_DIGITS_N + 1) * T_READ_CHAR

# Time to process one statement with the optimized algorithm.
# The algorithm reads 4 characters (the 3 statement chars + newline)
# and performs 1 character comparison.
time_per_statement = 4 * T_READ_CHAR + 1 * T_COMPARE_CHAR

# Total time to process all n statements in the loop.
total_time_statements = N_MAX * time_per_statement

# The final running time is the sum of reading 'n' and processing all statements.
# The time to print the final result is negligible (in nanoseconds).
total_time_ms = time_read_n + total_time_statements

# Step 4: Print the final calculation, showing each number involved.
print("Worst-case time estimation formula:")
print("T = (Time to read n) + (Number of statements * Time per statement)")
print("T = ((Digits in n + 1) * T_READ_CHAR) + (n * (4 * T_READ_CHAR + T_COMPARE_CHAR))")
print("\nSubstituting the values for the worst-case scenario (n=90):")
print(f"T = (({NUM_DIGITS_N} + 1) * {T_READ_CHAR}) + ({N_MAX} * (4 * {T_READ_CHAR} + {T_COMPARE_CHAR}))")
print(f"T = ({time_read_n}) + ({N_MAX} * ({time_per_statement}))")
print(f"T = {time_read_n} + {total_time_statements}")
print(f"T = {int(total_time_ms)} ms")