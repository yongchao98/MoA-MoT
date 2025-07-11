import math

def solve():
    """
    Calculates the minimal and maximal numbers of grid's cells that a
    circumference of radius 500 can cross.
    """
    R = 500

    # Step 1 & 2: Calculate the maximal number of cells
    # This occurs for a generic placement, e.g., center at (0.5, 0.5).
    # The formula is C_max = 8R.
    max_cells = 8 * R

    # Step 3, 4, 5: Calculate the minimal number of cells
    # This occurs for a special placement, e.g., center at (epsilon, epsilon).
    # The formula is C_min = 8R - 2(S - 4), where S = r_2(R^2).

    # First, calculate S = r_2(R^2).
    # R = 500 = 2^2 * 5^3
    # R^2 = 2^4 * 5^6
    # The prime factorization of R^2 has one prime of the form 4k+1, which is 5.
    # The exponent of 5 is 6.
    # r_2(n) = 4 * (a_1 + 1) * (a_2 + 1) * ... for primes p_i = 4k+1 with exponents a_i.
    # Here, we only have one such prime, 5, with exponent 6.
    s = 4 * (6 + 1)
    
    # The number of integer points on the circle x^2 + y^2 = 500^2 is S.
    num_points_on_circle = s
    
    # The number of these points not on the axes is S-4.
    # The four points on the axes are (R, 0), (-R, 0), (0, R), (0, -R).
    num_points_off_axes = num_points_on_circle - 4

    # Now, calculate the minimal number of cells.
    min_cells = 8 * R - 2 * num_points_off_axes

    # Step 7: Print the final answer in the specified format.
    print(f"The radius is R = {R}.")
    print(f"The maximal number of cells is calculated as 8 * R.")
    print(f"C_max = 8 * {R} = {max_cells}")
    print("\nTo find the minimum, we first calculate S = r_2(R^2), the number of integer points on x^2+y^2=R^2.")
    print(f"R^2 = 500^2 = 250000 = 2^4 * 5^6.")
    print(f"The number of points S = r_2(500^2) = 4 * (6 + 1) = {num_points_on_circle}.")
    print(f"The minimal number of cells is calculated as 8*R - 2*(S - 4).")
    print(f"C_min = 8 * {R} - 2 * ({num_points_on_circle} - 4) = {max_cells} - 2 * {num_points_off_axes} = {min_cells}")
    print(f"\nThe minimal and maximal numbers of cells are ({min_cells}, {max_cells}).")

solve()
<<< (3952, 4000) >>>