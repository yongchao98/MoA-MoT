import math
import sys

def gcd_list(numbers):
    """Computes the greatest common divisor of a list of numbers."""
    if not numbers:
        return 0
    # The GCD of a single number is the number itself.
    if len(numbers) == 1:
        return numbers[0]
    
    result = numbers[0]
    for i in range(1, len(numbers)):
        result = math.gcd(result, numbers[i])
    return result

def restore_sequence(s_str):
    """
    Restores a sequence A from a given set S.

    Args:
        s_str: A string containing space-separated integers representing the set S.
    
    Prints:
        The restored sequence A, or -1 if no such sequence exists.
    """
    try:
        # Step 1: Parse the input string into a sorted list of unique integers.
        s_set = set(map(int, s_str.split()))
        if not s_set:
            print("-1")
            return
        s_list = sorted(list(s_set))
    except (ValueError, IndexError):
        # Handle cases with invalid input format.
        print("-1")
        return

    # Step 2: Check the necessary and sufficient condition.
    # A solution exists if and only if the minimum element of S
    # is equal to the GCD of all elements in S.
    min_s = s_list[0]
    gcd_s = gcd_list(s_list)

    if min_s != gcd_s:
        # If the condition fails, no solution exists.
        print("-1")
        return

    # Step 3: If the condition holds, construct a valid sequence A.
    # A proven construction is to alternate the minimum element `g`
    # with the other elements of the set.
    if len(s_list) == 1:
        # If S has only one element, A is just that element.
        A = s_list
    else:
        A = []
        g = s_list[0]
        # Create the sequence [g, s2, g, s3, ...]
        for i in range(1, len(s_list)):
            A.append(g)
            A.append(s_list[i])
    
    # Step 4: Print the constructed sequence.
    # The problem asks to output the numbers in the final sequence.
    print(*A)

# Example usage with a sample input string.
# To use with user input, you would replace this with:
# input_s = sys.stdin.readline()
input_s = "6 12 18"
restore_sequence(input_s)