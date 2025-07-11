import sys

def solve():
    """
    This function calculates the number of non-admissible integers k
    for given positive integers a and b.
    """
    # We choose specific values for a and b as they are not provided.
    # The logic holds for any pair of positive integers.
    a = 4
    b = 7
    
    # Let's verify that a and b are positive integers
    if not (isinstance(a, int) and isinstance(b, int) and a > 0 and b > 0):
        print("Error: a and b must be positive integers.", file=sys.stderr)
        return

    # Case 1: a=1 or b=1
    if min(a, b) == 1:
        # The only admissible k is ab. Non-admissible are {0, ..., ab-1}.
        # The number of non-admissible integers is ab.
        result = a * b
        print(f"{a} * {b} = {result}")
    # Case 2: a>1 and b>1
    else:
        # The non-admissible integers are {1, ..., max(a,b)-1} and {ab-1}.
        # The number of non-admissible integers is max(a,b).
        result = max(a, b)
        print(f"max({a}, {b}) = {result}")

solve()