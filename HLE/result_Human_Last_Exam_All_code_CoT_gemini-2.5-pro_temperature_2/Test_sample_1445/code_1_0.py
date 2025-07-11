import math

def solve():
    """
    Calculates the minimum number of operations n needed to transform any given
    100-digit binary sequence into any other.

    The complexity of a sequence can be measured by its number of alternations
    (e.g., 0 to 1 or 1 to 0).
    - A sequence of 100 identical digits (e.g., "00...0") has 1 block and 0 alternations.
    - A sequence of 100 alternating digits (e.g., "0101...01") has 100 blocks and 99 alternations.

    An operation (inserting/removing a block of identical digits) can change the
    number of alternations by at most 2.
    - Insertion: d -> d d' d creates 2 new alternations.
    - Removal: d d' d -> d d destroys 2 alternations.

    The worst-case transformation is from a sequence with minimum complexity (0 alternations)
    to one with maximum complexity (99 alternations), or vice-versa.
    The number of operations is the total change required divided by the maximum change per operation.
    """

    # Total number of digits in the sequence
    num_digits = 100

    # Maximum number of alternations is num_digits - 1
    max_alternations = num_digits - 1

    # Maximum change in alternations per operation
    change_per_operation = 2

    # Calculate the minimum number of operations for the worst-case scenario
    n = math.ceil(max_alternations / change_per_operation)
    
    print("The maximum number of alternations to create or remove is 99.")
    print("Each operation can create or remove at most 2 alternations.")
    print("So, the number of operations is ceil(99 / 2).")
    print("99 / 2 = 49.5")
    print(f"ceil(49.5) = {int(n)}")
    
solve()