import math

# Plan:
# 1. Define the performance costs in milliseconds (ms) from the problem description.
# 2. Determine the maximum number of statements (n) that fit within the 366-character limit.
# 3. Calculate the time to read n in the worst case (i.e., for the maximum n).
# 4. Calculate the time to process a single statement in the worst case.
# 5. Sum the times to get the total estimated running time.

# Step 1: Define performance costs in milliseconds (ms)
T_READ_CHAR = 15
T_RECOGNIZE_DIGIT = 110
T_COMPARE_CHAR = 150
# Integer and print operations are negligible (nanoseconds), so their cost is 0.

# Step 2: Determine the maximum number of statements (n)
# A program's total length = len(str(n)) + 1 (for newline) + n * (3 chars/statement + 1 newline).
# Equation: len(str(n)) + 1 + 4*n <= 366.
# Let's assume n has 2 digits (e.g., n=90).
# Equation becomes: 2 + 1 + 4*n <= 366  =>  4n <= 363  => n <= 90.75.
# So, the maximum integer n is 90.
# Verification for n=90: 2 (for "90") + 1 (newline) + 90 * 4 (statements) = 3 + 360 = 363. This is <= 366.
max_n = 90
n_digits = len(str(max_n))

# Step 3: Calculate time to read n in the worst case (n=90)
# An efficient C loop to read n: while ((c = getchar()) != '\n') { ... }
# For each digit, we read, compare to '\n', and recognize the digit.
# For the final '\n', we read and compare.
time_per_digit = T_READ_CHAR + T_COMPARE_CHAR + T_RECOGNIZE_DIGIT
time_for_final_newline = T_READ_CHAR + T_COMPARE_CHAR
time_read_n = (n_digits * time_per_digit) + time_for_final_newline

# Step 4: Calculate time for one statement in the worst case
# An optimized interpreter checks the first character. The worst case is when it's 'X',
# requiring a second character read and more comparisons.
# Logic: if (c1=='+') {} else if (c1=='-') {} else { // 'X' case }
# Cost for "X++":
# 1. Read 'X': T_READ_CHAR
# 2. Compare 'X' vs '+': T_COMPARE_CHAR (fails)
# 3. Compare 'X' vs '-': T_COMPARE_CHAR (fails)
# 4. Read '+': T_READ_CHAR
# 5. Compare '+' vs '+': T_COMPARE_CHAR (succeeds)
# 6. Consume rest of line (one more '+' and one '\n'): 2 * T_READ_CHAR
time_statement_worst_case = (T_READ_CHAR + T_COMPARE_CHAR + T_COMPARE_CHAR +
                             T_READ_CHAR + T_COMPARE_CHAR +
                             T_READ_CHAR + T_READ_CHAR)

# Step 5: Calculate total time
total_time_statements = max_n * time_statement_worst_case
total_time = time_read_n + total_time_statements

# --- Final Calculation Output ---
print("Worst-Case Running Time Estimation")
print("="*40)

print(f"1. Time to read n (for worst-case n={max_n}):")
print(f"   ({n_digits} digits * ({T_READ_CHAR} + {T_COMPARE_CHAR} + {T_RECOGNIZE_DIGIT}) ms/digit) + ({T_READ_CHAR} + {T_COMPARE_CHAR}) ms for newline")
print(f"   = ({n_digits} * {time_per_digit}) + {time_for_final_newline} = {time_read_n} ms")
print("")

print(f"2. Time per statement (worst-case, e.g., 'X++'):")
print(f"   {T_READ_CHAR}(read) + {T_COMPARE_CHAR}(cmp) + {T_COMPARE_CHAR}(cmp) + {T_READ_CHAR}(read) + {T_COMPARE_CHAR}(cmp) + {T_READ_CHAR}(read) + {T_READ_CHAR}(read)")
print(f"   = {time_statement_worst_case} ms")
print("")

print(f"3. Total Time Calculation:")
print(f"   Total Time = Time(read n) + (n * Time(per statement))")
print(f"   Total Time = {time_read_n} + ({max_n} * {time_statement_worst_case})")
print(f"   Total Time = {time_read_n} + {total_time_statements}")
print(f"   Total Time = {int(total_time)} ms")
print("")

print("Final Answer:")
print(f"{int(total_time)} ms")