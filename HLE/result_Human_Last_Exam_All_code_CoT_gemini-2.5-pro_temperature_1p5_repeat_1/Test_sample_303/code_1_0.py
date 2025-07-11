import sys

# In the recipe, the function calls itself without a condition to stop.
# This would cause an infinite loop. We'll simulate this with a counter
# to prevent the script from actually crashing.
# Python's default recursion limit is usually 1000.
sys.setrecursionlimit(1500)
RECURSION_COUNTER = 0

def buggy_caramel_sauce():
    """Simulates the flawed recipe causing infinite recursion."""
    global RECURSION_COUNTER
    RECURSION_COUNTER += 1
    # This represents the "Serve with caramel sauce." line inside its own recipe
    buggy_caramel_sauce()

def fixed_caramel_sauce():
    """
    Simulates the recipe after applying the fix from option F.
    The self-call is removed, so the function simply terminates.
    """
    print("Making a single, finite batch of caramel sauce... Done.")

print("--- Simulating the Buggy Recipe ---")
print("Calling the 'Caramel Sauce' recipe as written leads to infinite recursion.")
try:
    buggy_caramel_sauce()
except RecursionError:
    print(f"Caught RecursionError after {RECURSION_COUNTER} calls. This is an infinite loop!")

print("\n--- Applying the Fix (Option F) ---")
print("Removing the recursive call fixes the infinite loop:")
fixed_caramel_sauce()

print("\n--- Correctly Calculating Fibonacci (The Recipe's Goal) ---")
memo = {}
def fibonacci(n):
    """A correct, efficient recursive Fibonacci function with a base case."""
    if n in memo:
        return memo[n]
    if n <= 1:
        return n
    # The actual recursive formula: F(n) = F(n-1) + F(n-2)
    memo[n] = fibonacci(n-1) + fibonacci(n-2)
    return memo[n]

target_num = 12
result = fibonacci(target_num)
prev_1 = fibonacci(target_num - 1)
prev_2 = fibonacci(target_num - 2)

print(f"To achieve the recipe's goal, a proper recursive function is needed.")
print(f"The equation for Fibonacci number {target_num} is F({target_num}) = F({target_num-1}) + F({target_num-2}).")
print("Here are the numbers in that final equation:")
print(f"{result} = {prev_1} + {prev_2}")