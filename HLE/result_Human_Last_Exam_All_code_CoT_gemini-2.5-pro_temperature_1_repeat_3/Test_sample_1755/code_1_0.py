import io
import sys

def solve_equilibrium():
    """
    Reads force vectors from stdin and determines if they are in equilibrium.
    """
    # This simulates the input from the problem description for demonstration.
    # In a real scenario, this would read from the actual standard input.
    test_input = """3
3 -1 7
-5 2 -4
2 -1 -3
"""
    sys.stdin = io.StringIO(test_input)

    try:
        # Read the number of vectors
        n = int(sys.stdin.readline())
    except (ValueError, IndexError):
        print("Invalid input: Could not read the number of vectors.")
        return

    # Store the components of each vector to display the full equation
    x_components = []
    y_components = []
    z_components = []

    # Read n vectors
    for _ in range(n):
        try:
            line = sys.stdin.readline()
            xi, yi, zi = map(int, line.split())
            x_components.append(xi)
            y_components.append(yi)
            z_components.append(zi)
        except (ValueError, IndexError):
            # Handle malformed lines if necessary
            continue

    # Calculate the sum of forces for each axis
    sum_x = sum(x_components)
    sum_y = sum(y_components)
    sum_z = sum(z_components)

    # Helper function to create a readable equation string like "3 + (-5) + 2 = 0"
    def format_equation(components, total):
        # Enclose negative numbers in parentheses for clarity
        parts = [f"({c})" if c < 0 else str(c) for c in components]
        return " + ".join(parts) + f" = {total}"

    # Print the final equations for each axis
    print(f"x-axis forces: {format_equation(x_components, sum_x)}")
    print(f"y-axis forces: {format_equation(y_components, sum_y)}")
    print(f"z-axis forces: {format_equation(z_components, sum_z)}")

    # A body is in equilibrium if the net force along each axis is zero
    if sum_x == 0 and sum_y == 0 and sum_z == 0:
        print("YES")
    else:
        print("NO")

# Execute the function
solve_equilibrium()