import math

# This script calculates the answer to Q4 and formats the final output string.

# --- Answers to Q1, Q2, Q3 ---
# Q1: 3 lines have compiler errors.
q1_answer = 3
# Q2: Minimal memory is 30D (12D for n, 12D for m, 6D for a).
q2_answer = 30
# Q3: Yes, the problem can be solved without branching.
q3_answer = "Y"

# --- Calculation for Q4 ---
n = 4 * 10**9
m = 4 * 10**9
a = 1

# Calculate the number of flagstones for each dimension using ceiling division.
# Python's math.ceil works, but in integer arithmetic it's (numerator + denominator - 1) // denominator
flagstones_n = math.ceil(n / a)
flagstones_m = math.ceil(m / a)

# Calculate the total number of flagstones.
q4_answer = flagstones_n * flagstones_m

# Print the equation for Q4 as requested.
# Note: The numbers are converted to integers to avoid scientific notation in the output.
print(f"Final equation for Q4: {int(flagstones_n)} * {int(flagstones_m)} = {int(q4_answer)}")

# --- Final Answer Formatting ---
# The final answer is presented in the format x:y:z:t
# where x,y,z,t are the answers to Q1, Q2, Q3, Q4 respectively.
final_answer_string = f"{q1_answer}:{q2_answer}:{q3_answer}:{int(q4_answer)}"
print(f"\n{final_answer_string}")