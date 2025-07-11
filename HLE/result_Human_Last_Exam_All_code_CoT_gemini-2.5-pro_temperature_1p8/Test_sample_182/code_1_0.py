# Constants from the problem description, converted to milliseconds.
TIME_READ_CHAR_MS = 15
TIME_RECOGNIZE_CHAR_MS = 110
MAX_PROGRAM_CHARS = 366

# --- Step 1: Find the maximum number of statements (n) ---
# The total characters in a program are (number of digits in n) + (n * 3 characters per statement).
# We must find the largest n such that `len(str(n)) + 3 * n <= 366`.
# By testing values, we find n=121 is the maximum:
# 3 digits for "121" + 3 * 121 statements = 3 + 363 = 366 characters.
MAX_N = 121
NUM_DIGITS_MAX_N = 3

# --- Step 2: Calculate execution time for different components ---

# Time to get a single recognized character from the tape (read + recognize).
# This is the most fundamental combined operation.
TIME_GET_RECOGNIZED_CHAR_MS = TIME_READ_CHAR_MS + TIME_RECOGNIZE_CHAR_MS

# Time for the worst-case statement (e.g., "X--").
# 1. Read and recognize 'X'.
# 2. Read and recognize '-'.
# 3. Read the final '-' to advance the tape past the line (no recognition needed).
TIME_PER_WORST_CASE_STATEMENT_MS = TIME_GET_RECOGNIZED_CHAR_MS + TIME_GET_RECOGNIZED_CHAR_MS + TIME_READ_CHAR_MS

# --- Step 3: Calculate the total worst-case execution time ---

# Time to read the number 'n' from the first line. In the worst case, n is 121,
# which has 3 digits that must each be read and recognized.
TIME_TO_READ_N_MS = NUM_DIGITS_MAX_N * TIME_GET_RECOGNIZED_CHAR_MS

# Total time for executing all 'n' statements, assuming every statement is the worst-case type.
TOTAL_TIME_FOR_STATEMENTS_MS = MAX_N * TIME_PER_WORST_CASE_STATEMENT_MS

# The final total time is the sum of parsing 'n' and processing all statements.
TOTAL_WORST_CASE_TIME_MS = TIME_TO_READ_N_MS + TOTAL_TIME_FOR_STATEMENTS_MS

# --- Step 4: Print the detailed calculation ---
print("Worst-Case Running Time Estimation")
print("-" * 35)
print(f"Max statements (n): {MAX_N}")
print(f"Time for one recognized char: {TIME_READ_CHAR_MS} + {TIME_RECOGNIZE_CHAR_MS} = {TIME_GET_RECOGNIZED_CHAR_MS} ms")
print(f"Worst-case time per statement: {TIME_GET_RECOGNIZED_CHAR_MS} + {TIME_GET_RECOGNIZED_CHAR_MS} + {TIME_READ_CHAR_MS} = {TIME_PER_WORST_CASE_STATEMENT_MS} ms")
print("-" * 35)
print("Final Calculation:")
print(f"Time to read n ('{MAX_N}'): {NUM_DIGITS_MAX_N} digits * {TIME_GET_RECOGNIZED_CHAR_MS} ms/digit = {TIME_TO_READ_N_MS} ms")
print(f"Time for {MAX_N} statements: {MAX_N} * {TIME_PER_WORST_CASE_STATEMENT_MS} ms/stmt = {TOTAL_TIME_FOR_STATEMENTS_MS} ms")
print(f"Total Estimated Time = {TIME_TO_READ_N_MS} + {TOTAL_TIME_FOR_STATEMENTS_MS}")
print(f"Result: {round(TOTAL_WORST_CASE_TIME_MS)} ms")