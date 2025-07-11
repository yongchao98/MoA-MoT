import sys

def solve_equilibrium():
    """
    Reads n 3D vectors, calculates the sum of forces on each axis,
    prints the equations, and determines if the body is in equilibrium.
    """
    try:
        # Read the number of vectors from the first line of input.
        n_line = sys.stdin.readline()
        if not n_line:
            return
        n = int(n_line.strip())

        # Store coordinates to print the final equations.
        x_coords, y_coords, z_coords = [], [], []
        
        # Loop n times to read all vectors.
        for _ in range(n):
            line = sys.stdin.readline()
            if not line:
                break # Stop if input ends prematurely.
            
            # Read the three integer coordinates from the line.
            xi, yi, zi = map(int, line.strip().split())
            
            x_coords.append(xi)
            y_coords.append(yi)
            z_coords.append(zi)

        # Calculate the sum of forces for each axis.
        sum_x = sum(x_coords)
        sum_y = sum(y_coords)
        sum_z = sum(z_coords)

        # Format the equations for each axis, wrapping numbers in parentheses for clarity.
        x_equation = ' + '.join(f'({i})' for i in x_coords) + f' = {sum_x}'
        y_equation = ' + '.join(f'({i})' for i in y_coords) + f' = {sum_y}'
        z_equation = ' + '.join(f'({i})' for i in z_coords) + f' = {sum_z}'
        
        # Print the final equations.
        print(x_equation)
        print(y_equation)
        print(z_equation)

        # Check if the sum of forces on all axes is zero.
        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            print("NO")

    except (ValueError, IndexError):
        # Handle cases with malformed input.
        print("NO")

# Execute the solution function.
solve_equilibrium()