import sys

def check_equilibrium():
    """
    Reads force vectors from standard input, calculates the resultant force,
    and determines if the body is in equilibrium.
    """
    try:
        # Read the number of force vectors.
        num_vectors = int(sys.stdin.readline())
    except (ValueError, IndexError):
        # Handle empty or invalid input for the number of vectors.
        # A robust program should not crash on bad input.
        print("Invalid input: Please provide a valid number of vectors.")
        return

    # Initialize the sum of force components for each axis.
    sum_x, sum_y, sum_z = 0, 0, 0

    # Read each vector and accumulate the components.
    for _ in range(num_vectors):
        try:
            line = sys.stdin.readline()
            # If input ends prematurely, stop.
            if not line:
                break
            # Split the line into components and add to the sums.
            xi, yi, zi = map(int, line.split())
            sum_x += xi
            sum_y += yi
            sum_z += zi
        except (ValueError, IndexError):
            # Skip malformed lines.
            continue
            
    # The final equation for equilibrium is that the sum of force vectors is zero.
    # The numbers in this equation are the components of the sum vector.
    # As per instructions, printing these numbers.
    print(f"Resultant force components: Fx={sum_x}, Fy={sum_y}, Fz={sum_z}")
    
    # Check for equilibrium condition: the resultant force vector must be (0, 0, 0).
    if sum_x == 0 and sum_y == 0 and sum_z == 0:
        print("YES")
    else:
        print("NO")

if __name__ == "__main__":
    check_equilibrium()