from re import *

# The 32-character string to be placed in place of @ is:
# *(sorted({int(i)for i in findall(

# The complete line of code would be:
# print(*(sorted({int(i)for i in findall(r'\d+', input())}))or ["NO"])

# Here is the code to demonstrate the solution.
# I will simulate the one-liner by building the expression string and executing it.

# Step 1: Define the missing part
replacement_string = "*(sorted({int(i)for i in findall("

# Step 2: Get user input for the demonstration
user_input = input()

# Step 3: Construct the full command as a single string
# We need to pass the user_input securely into the eval context
# to prevent injection issues, although for this specific problem it's less of a risk.
# The `namespace` dictionary provides a safe context for evaluation.
code_to_run = f'print({replacement_string}r\\\'\\d+\\\', a)})}}))or ["NO"])'
# A slightly cleaner way to write this without complex escaping for the print statement
code_to_run = 'print(*(sorted({int(i) for i in findall(r\'\\d+\', user_input)})) or ["NO"])'


# Step 4: Execute the code using eval
# Using exec as it's a statement. globals() and a local dict are passed for context.
# It's generally unsafe to use exec/eval on user input, but here we control the structure
# and the input is only used within findall.
exec(code_to_run, {'findall': findall, 'user_input': user_input})
