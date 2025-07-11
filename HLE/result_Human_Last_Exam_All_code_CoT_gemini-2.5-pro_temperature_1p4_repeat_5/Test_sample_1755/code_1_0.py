import sys

def solve_and_explain():
    """
    This function reads vector data from standard input, calculates the resultant force,
    prints the breakdown of the sum for each component, and determines if a body
    is in equilibrium.
    """
    try:
        # Read the number of vectors, n
        n_line = sys.stdin.readline()
        if not n_line.strip():
            # Handle empty input
            return
        n = int(n_line.strip())

        # Ensure n is within the specified constraints
        if not (1 <= n <= 100):
            return

        x_coords = []
        y_coords = []
        z_coords = []

        # Read n vectors from the following lines
        for _ in range(n):
            line = sys.stdin.readline()
            if not line.strip():
                # Stop if input ends prematurely
                break
            xi, yi, zi = map(int, line.strip().split())
            x_coords.append(xi)
            y_coords.append(yi)
            z_coords.append(zi)
            
        # If not enough vectors were read, it's an invalid input.
        if len(x_coords) != n:
            return

        # Calculate the sum for each coordinate
        sum_x = sum(x_coords)
        sum_y = sum(y_coords)
        sum_z = sum(z_coords)

        # Print the equation for the sum of x-components
        # The map(str, ...) converts each integer number to a string for joining
        print(" + ".join(map(str, x_coords)) + f" = {sum_x}")

        # Print the equation for the sum of y-components
        print(" + ".join(map(str, y_coords)) + f" = {sum_y}")

        # Print the equation for the sum of z-components
        print(" + ".join(map(str, z_coords)) + f" = {sum_z}")

        # Check for equilibrium: the sum of vectors is the zero vector
        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            print("NO")

    except (ValueError, IndexError):
        # Handles cases with non-integer input or improperly formatted lines
        return

solve_and_explain()