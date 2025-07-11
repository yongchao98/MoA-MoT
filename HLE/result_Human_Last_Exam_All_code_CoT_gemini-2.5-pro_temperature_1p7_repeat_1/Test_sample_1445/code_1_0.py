def solve():
    """
    Calculates the minimum number of operations n needed to transform any given
    100-digit sequence of 0s and 1s into any other.

    The plan is as follows:
    1.  Identify the worst-case transformation. This is likely from a sequence of
        minimal complexity (e.g., all identical digits, 1 run) to one of
        maximal complexity (e.g., fully alternating digits, 100 runs).
        Let's choose the initial sequence Si = '0' * 100 and the target
        sequence St = '01' * 50.

    2.  Calculate the cost of a plausible transformation path: Si -> Smid -> St.
        a.  Transform Si to an intermediate sequence Smid = '0' * 50 + '1' * 50.
            - Step 1: '0'*100 -> '0'*50. This is one 'delete' operation.
            - Step 2: '0'*50 -> '0'*50 + '1'*50. This is one 'insert' operation.
            - Cost for this stage is 2.

        b.  Transform Smid to St. This involves interleaving the 50 '0's and 50 '1's.
            - Smid has 2 runs. St has 100 runs. We need to create 98 new run boundaries.
            - A "splitting insertion" (e.g., '000' -> '010') adds 2 runs in one operation.
            - To go from 2 runs to 100 runs, we need (100 - 2) / 2 = 49 splitting insertions.
            - These 49 insertions create a sequence with the correct run structure
              but the wrong length (100 + 49 = 149 digits). We need to remove the
              49 extra digits. This requires 49 deletion operations.
            - Cost for this stage is 49 insertions + 49 deletions = 98.

    3.  Sum the costs.
        - Total cost = (cost of Si -> Smid) + (cost of Smid -> St)
        - Total cost = 2 + 98 = 100.
    """
    # Cost to transform '0'*100 into '0'*50 + '1'*50
    cost_to_mid = 1 + 1  # 1 deletion, 1 insertion
    
    # Cost to interleave '0'*50 + '1'*50 into '01'*50
    # This involves creating 98 new alternating runs from 2 runs.
    # Each "split" operation creates 2 new runs. Number of splits = (100 - 2) / 2 = 49.
    # Each split is one insertion, adding a digit. 49 insertions.
    # We now have 100+49=149 digits, need to remove 49. 49 deletions.
    cost_to_interleave = 49 + 49

    # The maximum number of operations is the sum of these costs.
    n = cost_to_mid + cost_to_interleave
    
    print("The minimum number of operations n needed is the sum of the steps:")
    print(f"1. Cost to create a sequence with the right number of 0s and 1s ('0'*50 + '1'*50): {cost_to_mid}")
    print(f"2. Cost to interleave the digits to form the target run-sequence: {cost_to_interleave}")
    print(f"Final equation: {cost_to_mid} + {cost_to_interleave} = {n}")
    
solve()