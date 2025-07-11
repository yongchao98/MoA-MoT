import math

def get_divisors_greater_than_one(n):
    """Returns a set of all divisors of n that are greater than 1."""
    if n <= 1:
        return set()
    divs = set()
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    divs.add(n)
    return divs

# From Sharkovsky's Theorem and the problem's conditions, the set of minimal
# periods that do not exist for the function f is {3, 5, 7, 9, 11}.
non_existent_minimal_periods = {3, 5, 7, 9, 11}

# The set S contains integers k for which there is no point of order k.
# This means that for any k in S, all of its divisors d > 1 must be
# in the set of non_existent_minimal_periods.
S = set()

# We can find all elements of S by checking integers up to a reasonable bound.
# Any k with a divisor not in non_existent_minimal_periods will be excluded.
for k in range(1, 100):
    divisors_k = get_divisors_greater_than_one(k)
    if divisors_k.issubset(non_existent_minimal_periods):
        S.add(k)

# The prompt asks to output each number in the "final equation",
# which we interpret as the elements of the set S.
print("The elements of the set S are:")
# We print the sorted list for clarity.
final_set_S = sorted(list(S))
for number in final_set_S:
    print(number)

print("\nThe cardinality of S is:")
print(len(final_set_S))