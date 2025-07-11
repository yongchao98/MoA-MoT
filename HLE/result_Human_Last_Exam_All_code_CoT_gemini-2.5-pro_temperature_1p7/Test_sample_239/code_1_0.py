import math

# Answer to the four questions about the XVM program
# 1. How many lines of code have compiler errors?
#    - Line 3: `unsigned long long` is not a valid type in XVM.
#    - Line 4: `scanf` uses `%d` (for type `digit`) with pointers to `unsigned long long`. The format specifier is incorrect. It should be `%n` for `unsigned long`.
#    - Line 9: `printf` uses `%d` (for type `digit`) to print a large number. The format specifier is incorrect. It should be `%n`.
#    So, there are 3 lines with errors.

# 2. Rewrite it to use the least amount of memory and the smallest number of statements with branching. What is the total memory (in D) used for variables?
#    - To avoid branching (if statements), we can calculate the ceiling of a division `x/y` using integer arithmetic: `(x + y - 1) / y`.
#    - We need variables for n, m, and a.
#    - n and m can be up to 4*10^9, which requires the `unsigned long` type (12D).
#    - a can be up to 10^6, which also fits in `unsigned long`. Using the same type for all variables is best for calculations.
#    - Minimal variables needed: n, m, a. Total memory: 12D (for n) + 12D (for m) + 12D (for a) = 36D.

# 3. Can this problem be solved without branching statements?
#    - Yes. Using the integer arithmetic `(x + y - 1) / y` to calculate the ceiling avoids `if` statements.

# 4. What is the expected output for m = n = 4*10^9 and a = 1.
#    - na = ceil(4*10^9 / 1) = 4*10^9
#    - ma = ceil(4*10^9 / 1) = 4*10^9
#    - Total = 4*10^9 * 4*10^9 = 16 * 10^18 = 16000000000000000000.

# The Python code below solves the Theatre Square problem itself.
def solve():
    try:
        # Read n, m, and a from a single line of input
        n, m, a = map(int, input().split())

        # Calculate the number of flagstones for the 'n' side using ceiling division
        # (n + a - 1) // a is a way to do ceiling division with integers
        flagstones_n = (n + a - 1) // a

        # Calculate the number of flagstones for the 'm' side
        flagstones_m = (m + a - 1) // a

        # Total flagstones is the product of the two
        total_flagstones = flagstones_n * flagstones_m

        # Print the final equation
        print(f"{flagstones_n} * {flagstones_m} = {total_flagstones}")

    except (ValueError, IndexError):
        print("Please provide three positive integers separated by spaces (e.g., 6 7 4)")

# The user is expected to run this function and provide input.
# For example, if the input is "4000000000 4000000000 1", the code will print the required result.
# To demonstrate, we can call it here for testing, but typically the user would execute it.
# To match the interactive format, the function is defined but not called directly in the final script.
# If you were to run this file, you would call solve()
# solve()
# The final formatted answer is provided after the code block.
