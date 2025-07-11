import sys

# This script generates the Tate-style program-expression graph for the given loop.

# 1. Define the components of the μ-expression based on the C-like code.
# The code is:
# int a = 1;
# for (int i = 1; i < n; i++) { a *= i; }

# The bound variable, as requested by the prompt.
bound_variable = "a"

# The initial value of the variable 'a' before the loop.
initial_value = 1

# The update expression inside the loop. The bound variable 'a' represents
# the value from the previous iteration. 'i' is the loop counter.
update_expression = "a * i"

# 2. Construct the final expression string in the format μ x.(initial, update).
# We explicitly include the number `1` from the initial value in the final string.
final_expression = f"μ {bound_variable}.({initial_value}, {update_expression})"

# 3. Print the final expression.
print(final_expression)

# Adding the answer in the specified format for clarity, as requested.
# This part is technically metadata and the primary output is from the print statement above.
sys.stdout.flush() # ensure the print is visible before the final answer tag.
# The <<<...>>> format is for the final, concise answer extraction.
final_answer_for_extraction = final_expression
# print(f"\n<<<{final_answer_for_extraction}>>>")