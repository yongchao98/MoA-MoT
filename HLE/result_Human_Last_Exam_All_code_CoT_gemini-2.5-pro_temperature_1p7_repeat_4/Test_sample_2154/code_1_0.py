import math

def get_ur_n(n):
    """
    Calculates the minimal order u_r(n) based on whether n is even or odd.
    """
    if n % 2 != 0:
        # For odd n, the order is n-1
        return n - 1
    else:
        # For even n, the order is 2 * floor(n/4) due to symmetry
        return 2 * math.floor(n / 4)

def solve_problem():
    """
    Calculates and prints the sequence of u_r(n) for n from 3 to 12.
    """
    print("The values for {u_r(3), u_r(4), ..., u_r(12)} are:")
    for n in range(3, 13):
        ur_val = get_ur_n(n)
        print(f"u_r({n}) = {ur_val}")

if __name__ == "__main__":
    solve_problem()