import math

def count_maximal_cool_strings(n):
    """
    Calculates and explains the number of 'cool strings' of maximal length for n symbols.

    Args:
        n (int): The number of distinct symbols. Must be a non-negative integer.
    """
    if not isinstance(n, int) or n < 0:
        print("Error: The number of symbols (n) must be a non-negative integer.")
        return

    print(f"For n = {n} different symbols, we want to find the number of 'cool strings' of maximal length, which is 3n = {3*n}.")
    print("\nBased on the analysis, a cool string of maximal length is constructed by:")
    print("1. Taking a permutation (an ordered arrangement) of the n symbols.")
    print("2. Creating the string by concatenating each symbol three times according to that permutation.")
    print("\nFor example, if n=2 with symbols A and B, the possible permutations are (A, B) and (B, A).")
    print("This corresponds to the cool strings 'AAABBB' and 'BBBAAA'.")
    print("\nThe total number of such strings is the number of permutations of n symbols, which is n factorial (n!).")

    if n == 0:
        result = 1
        equation = f"{n}! = 1"
    else:
        result = math.factorial(n)
        # Showing the expanded form of the factorial calculation
        equation_parts = [str(i) for i in range(1, n + 1)]
        equation = f"{n}! = {' * '.join(equation_parts)}"

    print(f"\nThe calculation for n = {n} is:")
    if n > 0:
      print(f"{equation} = {result}")
    else:
      print(equation)


# --- Main execution ---
# The user can change this value to find the answer for a different number of symbols.
n_value = 5
count_maximal_cool_strings(n_value)
# The final result for n=5 is 120. The numeric part of the answer is extracted below.
result_for_n_5 = math.factorial(n_value)
print(f"\n<<<For n={n_value}, the answer is {result_for_n_5}>>>")