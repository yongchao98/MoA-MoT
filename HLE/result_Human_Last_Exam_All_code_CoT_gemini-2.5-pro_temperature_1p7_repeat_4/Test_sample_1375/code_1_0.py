import math

def double_factorial(n):
    """Computes the double factorial of n."""
    if n < 0:
        return 0
    if n == 0:
        return 1
    res = 1
    i = n
    while i >= 2:
        res *= i
        i -= 2
    return res

def solve():
    """
    Calculates the average area of a projection of the n-dimensional cross-polytope
    for odd dimension n = 2k + 1.
    The user is prompted to enter n.
    """
    try:
        n_str = input("Enter an odd dimension n (e.g., 3, 5, 7, ...): ")
        n = int(n_str)
        if n <= 0 or n % 2 == 0:
            print("Dimension n must be a positive odd integer.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    # k = (n - 1) // 2
    # The formula is (n-2)!! / (n-3)!!
    # Let's handle the n=1 case separately as (1-3)!! is not standard
    if n == 1:
        # A 1D cross-polytope is the interval [-1, 1]. A projection on a hyperplane (a point) has area 2 (the two points {-1,1}).
        # The average area after rotation does not make much sense, but using the formula gives ( -1)!! / (-2)!! = 1/0, not defined
        # However, a "projection" of a line segment onto a point is just the endpoints if the line is orthogonal to the "hyperplane" (point).
        # The concept of random rotation is also tricky. The limit of our formula for n->1 doesn't seem to work.
        # But if we consider a line segment in R^1, its only projection is onto {0}, which is just a point (area 0?).
        # However, for a 2D cross-polytope (a square rotated on the plane), the projection is a line segment, whose length (1-area) is averaged.
        # Average length is 2*sqrt(2)*2/pi for square of side sqrt(2). Our crosspolytope is rotated diamond, side 1. So 4/pi.
        # The formula only makes sense for n >= 3
        if n==1:
            print("The problem is ill-defined for n=1. The formula applies for n>=3.")
            print("For n = 1, k=0. The formula would be (-1)!!/(-2)!! which is not well-defined using standard double factorials.")
            print("Let's compute for n >= 3.")
            return

    numerator = double_factorial(n - 2)
    denominator = double_factorial(n - 3)
    
    if denominator == 0:
        print(f"Error: division by zero for n={n}")
        return

    result = numerator / denominator

    print(f"For dimension n = {n}:")
    print(f"The average area is calculated using the formula (n-2)!! / (n-3)!!")
    print(f"A_avg = {n-2}!! / {n-3}!! = {numerator} / {denominator} = {result}")

solve()