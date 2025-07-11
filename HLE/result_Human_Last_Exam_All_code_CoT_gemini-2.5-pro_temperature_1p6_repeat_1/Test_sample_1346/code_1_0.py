def calculate_and_print_a(n, p):
    """
    Calculates a(n) using the recurrence relation a_k = 4*a_{k-1} - a_{k-2}.
    It also prints the steps of the calculation as requested.
    """
    if n == 0:
        print(f"For p={p}, we calculate a(0) = 1")
        return 1
    if n == 1:
        print(f"For p={p}, we calculate a(0)=1, a(1)=3")
        return 3

    print(f"For p={p}, we need to calculate a({n}). We start with a(0)=1, a(1)=3.")
    a0, a1 = 1, 3
    print(f"a(0) = {a0}")
    print(f"a(1) = {a1}")
    for i in range(2, n + 1):
        an = 4 * a1 - a0
        print(f"a({i}) = 4 * a({i-1}) - a({i-2}) = 4 * {a1} - {a0} = {an}")
        a0, a1 = a1, an
    return a1

# Case 1: p = 50051
p1 = 50051
# The index N mod (p-1) is 5
n1 = 5
result1 = calculate_and_print_a(n1, p1)

print("-" * 20)

# Case 2: p = 50069
p2 = 50069
# The index N mod (p+1) is 3
n2 = 3
result2 = calculate_and_print_a(n2, p2)

print("-" * 20)
print("The calculated values are:")
print(f"{result1},{result2}")
