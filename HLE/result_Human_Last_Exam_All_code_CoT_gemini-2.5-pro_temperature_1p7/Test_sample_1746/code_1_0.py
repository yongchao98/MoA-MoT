# This script analyzes and demonstrates the correctness of statement (B) regarding
# Python's division semantics evolution.

# Statement B: The expression 1/2 * 2 equals 1/2 + 1/2 in all Python versions.

# --- Analysis ---
# This statement is true across Python's major versions, but the underlying
# numerical values are different due to the change in the `/` operator.

# 1. In Python 2, `/` with integers performs integer division.
#    - `1 / 2` evaluates to 0.
#    - The left side becomes `0 * 2`, which is 0.
#    - The right side becomes `0 + 0`, which is 0.
#    - The equality `0 == 0` is true.

# 2. In Python 3, `/` performs "true" or float division.
#    - `1 / 2` evaluates to 0.5.
#    - The left side becomes `0.5 * 2`, which is 1.0.
#    - The right side becomes `0.5 + 0.5`, which is 1.0.
#    - The equality `1.0 == 1.0` is true.

# The code below demonstrates these two cases.

print("--- Demonstrating in Python 3 (True Division) ---")
# The original numbers in the expression
num1 = 1
num2 = 2
val_1_div_2 = num1 / num2
left_side = val_1_div_2 * num2
right_side = val_1_div_2 + val_1_div_2
print(f"Expression: ({num1}/{num2} * {num2}) == ({num1}/{num2} + {num1}/{num2})")
print(f"Evaluated: ({val_1_div_2} * {num2}) == ({val_1_div_2} + {val_1_div_2})")
print(f"Final Values: {left_side} == {right_side}")
print(f"Is the statement true in Python 3? {left_side == right_side}")

print("\n--- Simulating Python 2 Behavior (Integer Division) ---")
# In Python 3, `//` performs the integer division that `/` did in Python 2.
val_1_div_2_sim = num1 // num2
left_side_sim = val_1_div_2_sim * num2
right_side_sim = val_1_div_2_sim + val_1_div_2_sim
print(f"Simulating Py2: ({num1}//{num2} * {num2}) == ({num1}//{num2} + {num1}//{num2})")
print(f"Evaluated: ({val_1_div_2_sim} * {num2}) == ({val_1_div_2_sim} + {val_1_div_2_sim})")
print(f"Final Values: {left_side_sim} == {right_side_sim}")
print(f"Is the statement true for Python 2's behavior? {left_side_sim == right_side_sim}")
