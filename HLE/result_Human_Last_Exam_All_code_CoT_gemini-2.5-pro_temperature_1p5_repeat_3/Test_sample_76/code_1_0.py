import math

def is_perfect_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False, -1
    if n == 0:
        return True, 0
    sqrt_n = int(math.sqrt(n))
    if sqrt_n * sqrt_n == n:
        return True, sqrt_n
    return False, -1

def solve_equation():
    """
    Finds the number of integers x for which x^3 - 16x^2 - 72x + 1056 is a perfect square.
    The search is narrowed down to the ranges [-8, 7] and [17, 62].
    """
    solutions = []
    # Combine the search ranges
    search_ranges = list(range(-8, 8)) + list(range(17, 63))

    for x in search_ranges:
        # Define the polynomial P(x)
        val = x**3 - 16*x**2 - 72*x + 1056
        is_sq, y = is_perfect_square(val)
        if is_sq:
            solutions.append(x)
            print(f"For x = {x}, the equation is:")
            print(f"{x}^3 - 16*({x})^2 - 72*({x}) + 1056 = {val} = {y}^2")
            print("-" * 20)

    print(f"The total number of integers x for which the quantity is a perfect square is {len(solutions)}.")
    return len(solutions)

if __name__ == '__main__':
    num_solutions = solve_equation()
    # The final answer is wrapped according to the required format.
    # The actual numerical answer is returned by the function for clarity,
    # but we print it here in the requested final format.
    # print(f"<<<{num_solutions}>>>") # This would be the final output line if run in a special environment
                                     # For this environment, let's just output the number.
    print(num_solutions)
