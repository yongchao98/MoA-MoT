import sys

def solve():
    """
    Reads force vectors, checks for equilibrium, and prints the process.
    """
    try:
        # Read the number of vectors, n.
        n_line = sys.stdin.readline()
        if not n_line.strip():
            print("NO")
            return
        n = int(n_line)
    except (ValueError, IndexError):
        # Handle cases where n is not a valid integer or the line is empty.
        print("NO")
        return

    if n == 0:
        print("YES")
        return

    sum_x, sum_y, sum_z = 0, 0, 0
    x_coords, y_coords, z_coords = [], [], []

    for _ in range(n):
        try:
            line = sys.stdin.readline()
            if not line.strip():
                # Stop if we run out of input lines prematurely.
                break
            parts = list(map(int, line.split()))
            # Ensure the line has exactly three components.
            if len(parts) != 3:
                continue
            
            xi, yi, zi = parts
            
            x_coords.append(xi)
            y_coords.append(yi)
            z_coords.append(zi)
            
            sum_x += xi
            sum_y += yi
            sum_z += zi
        except (ValueError, IndexError):
            # Skip malformed lines.
            continue
    
    # Format and print the summation equations.
    # The str() conversion is important to handle negative numbers correctly.
    print(' + '.join(map(str, x_coords)) + f" = {sum_x}")
    print(' + '.join(map(str, y_coords)) + f" = {sum_y}")
    print(' + '.join(map(str, z_coords)) + f" = {sum_z}")

    # Check for equilibrium and print the final result.
    if sum_x == 0 and sum_y == 0 and sum_z == 0:
        print("YES")
    else:
        print("NO")

solve()