import math

def gcd(a, b):
    """Computes the greatest common divisor of two positive integers."""
    while b:
        a, b = b, a % b
    return abs(a)

def gcd_list(numbers):
    """Computes the greatest common divisor of a list of integers."""
    if not numbers:
        return 0
    result = numbers[0]
    for i in range(1, len(numbers)):
        result = gcd(result, numbers[i])
    return result

def solve_and_print(s_elements_str):
    """
    Reads a string of numbers representing set S, solves, and prints the result.
    """
    s_elements = list(map(int, s_elements_str.split()))

    if not s_elements:
        print("-1")
        return

    # 1. Find the minimum element of S
    min_s = min(s_elements)

    # 2. Calculate the GCD of all elements in S
    gcd_s = gcd_list(s_elements)

    # 3. Check the necessary and sufficient condition
    if min_s != gcd_s:
        print("-1")
        return

    # 4. If the condition holds, construct and print a valid sequence A
    g = min_s
    other_elements = [s for s in s_elements if s != g]

    if not other_elements:
        # Case where S has only one element, e.g., S = {5}
        reconstructed_A = [g]
    else:
        # Interleave other elements with g: A = [s'_1, g, s'_2, ...]
        reconstructed_A = []
        # Start with the first element from the 'others' list
        reconstructed_A.append(other_elements[0])
        # Loop through the rest of the 'others' list, interleaving with g
        for i in range(1, len(other_elements)):
            reconstructed_A.append(g)
            reconstructed_A.append(other_elements[i])

    # Print the length of the reconstructed sequence A
    print(len(reconstructed_A))
    # Print the elements of the reconstructed sequence A
    print(*reconstructed_A)


# The user did not provide input, so we use a placeholder example
# to demonstrate the code's functionality.
# Example Input: A set S = {2, 4, 8}
# min(S) = 2, gcd(2, 4, 8) = 2. Condition holds.
# g = 2, other_elements = [4, 8]
# Reconstructed A = [4, 2, 8]
# Output:
# 3
# 4 2 8
placeholder_input = "2 4 8"
solve_and_print(placeholder_input)
