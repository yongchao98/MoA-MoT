# Plan: Calculate the worst-case running time for an optimized X++ interpreter.
# The interpreter reads char by char to find the operation ('+' or '-'),
# then skips the rest of the line.

# --- Performance Timings in milliseconds (ms) ---
T_READ_CHAR = 15           # Time to read a character from tape
T_RECOGNIZE_CHAR = 110     # Time to recognize a character
T_EOLN_CHECK = 15          # Time to check for end of line (assumed to be a read)

# --- Derived Timings ---
# Time to read and recognize a character
T_READ_AND_RECOGNIZE = T_READ_CHAR + T_RECOGNIZE_CHAR
# Time to read a character just to skip it
T_SKIP_CHAR = T_READ_CHAR

# --- Step 1: Calculate the maximum number of statements (n) ---
# Total characters on tape <= 366
# Program structure: [n_digits] [newline] [n statements of 3 chars] [n newlines]
# Total chars = len(str(n)) + 1 + n * (3 + 1)
# For a 2-digit n, len(str(n)) = 2.
# 2 + 1 + 4*n <= 366  =>  4n <= 363  => n <= 90.75
max_n = 90

# --- Step 2: Calculate time to read n in the worst case (n=90) ---
# Input on tape: "90\n"
# This requires reading and recognizing two digits.
# A read operation is preceded by an eoln() check.
time_to_read_digit = T_EOLN_CHECK + T_READ_AND_RECOGNIZE
num_digits_in_max_n = 2
time_to_read_n_digits = num_digits_in_max_n * time_to_read_digit

# After reading the digits, a final eoln() check returns true, terminating the read loop.
final_eoln_check_for_n = T_EOLN_CHECK
time_to_read_n = time_to_read_n_digits + final_eoln_check_for_n

# --- Step 3: Calculate time for the worst-case statement ("X++\n") ---
# The interpreter must perform two read-and-recognize operations.
# 1. Read 'X': eoln() check (false) + read/recognize 'X'
time_for_first_char = T_EOLN_CHECK + T_READ_AND_RECOGNIZE
# 2. Read '+': eoln() check (false) + read/recognize '+'
time_for_second_char = T_EOLN_CHECK + T_READ_AND_RECOGNIZE
#    Operation is now identified. The rest of the line is skipped.
# 3. Skip '+': eoln() check (false) + read (skip) '+'
time_to_skip_third_char = T_EOLN_CHECK + T_SKIP_CHAR
# 4. Final eoln() check returns true, terminating the statement read.
final_eoln_check_for_statement = T_EOLN_CHECK

time_per_statement_worst = time_for_first_char + time_for_second_char + time_to_skip_third_char + final_eoln_check_for_statement

# --- Step 4: Calculate total worst-case running time ---
# Total Time = (Time to read n) + n * (Time for one worst-case statement)
total_time = time_to_read_n + max_n * time_per_statement_worst

# --- Final Output ---
print("Worst-Case Running Time Estimation:")
print(f"1. Max number of statements (n): {max_n}")
print(f"2. Time to read n={max_n} (2 digits):")
print(f"   ({T_EOLN_CHECK} + {T_READ_AND_RECOGNIZE}) * {num_digits_in_max_n} digits + {final_eoln_check_for_n} ms (final EOLN check) = {time_to_read_n} ms")
print(f"3. Time for one worst-case statement (e.g., 'X++'):")
print(f"   ({T_EOLN_CHECK} + {T_READ_AND_RECOGNIZE}) for 'X' + ({T_EOLN_CHECK} + {T_READ_AND_RECOGNIZE}) for '+' + ({T_EOLN_CHECK} + {T_SKIP_CHAR}) to skip last char + {final_eoln_check_for_statement} ms (final EOLN check) = {time_per_statement_worst} ms")
print("\nFinal Calculation:")
print(f"Total Time = Time to read n + n * Time per statement")
print(f"Total Time = {time_to_read_n} ms + {max_n} * {time_per_statement_worst} ms = {int(total_time)} ms")
