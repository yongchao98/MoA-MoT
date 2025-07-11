import math

def solve_circle_grid_problem():
    """
    Calculates the minimal and maximal numbers of grid cells a circle can cross.
    """
    R = 500

    # The maximum number of cells crossed is 8R for a generic circle center.
    max_cells = 8 * R

    # To find the minimum, we need to count the number of non-trivial integer
    # solutions to a^2 + b^2 = R^2, which are Pythagorean triples.
    # L is the number of pairs (a, b) with 0 < a < b and a^2 + b^2 = R^2.
    L = 0
    # We only need to check for 'a' up to R/sqrt(2) to ensure a < b.
    limit = int(R / math.sqrt(2))
    
    found_pairs = []
    for a in range(1, limit + 1):
        b_squared = R**2 - a**2
        if b_squared > 0:
            b = math.isqrt(b_squared)
            if b * b == b_squared:
                # Found an integer solution (a, b). Since a <= limit, a < b is guaranteed
                # because a=b only if a = R/sqrt(2), which is not an integer.
                L += 1
                found_pairs.append((a, b))

    # The minimal number of cells is 8R - 4L.
    min_cells = 8 * R - 4 * L
    
    print(f"Radius R = {R}")
    print(f"Number of Pythagorean triples (a, b, R) with 0 < a < b is L = {L}.")
    print("The found pairs (a, b) are:", found_pairs)
    print("\nCalculating the minimal and maximal number of cells:")
    print(f"Maximal number = 8 * {R} = {max_cells}")
    print(f"Minimal number = 8 * {R} - 4 * {L} = {min_cells}")
    
    print(f"\nThe final answer is ({min_cells}, {max_cells}).")


solve_circle_grid_problem()
<<< (3988, 4000) >>>