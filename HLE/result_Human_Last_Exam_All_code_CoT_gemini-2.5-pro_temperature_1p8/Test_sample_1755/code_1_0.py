import sys

def solve():
    """
    Reads force vectors and determines if a body is in equilibrium.
    """
    try:
        lines = sys.stdin.readlines()
        if not lines:
            print("NO")
            return

        # The first line should contain n, but the buggy C code ignores it
        # and loops until EOF. We will process all vector lines provided.
        # Assuming the first line is n and we should skip it to match the problem spec.
        try:
            n = int(lines[0].strip())
            vector_lines = lines[1:1+n]
        except (ValueError, IndexError):
            # If the first line is not a valid n, treat all lines as vectors
            vector_lines = lines

        x_sum, y_sum, z_sum = 0, 0, 0
        
        for line in vector_lines:
            parts = line.split()
            if len(parts) == 3:
                try:
                    xi = int(parts[0])
                    yi = int(parts[1])
                    zi = int(parts[2])
                    x_sum += xi
                    y_sum += yi
                    z_sum += zi
                except ValueError:
                    # Ignore malformed lines
                    continue
        
        # The correct logic for the C program is to sum the force vectors.
        # The body is in equilibrium if the sum of forces on each axis is zero.
        # So we check if x_sum, y_sum, and z_sum are all zero.
        # Final equation: x_sum + y_sum + z_sum = 0
        # More precisely, we need x_sum=0, y_sum=0, and z_sum=0.
        
        print(f"Sum of x components: {x_sum}")
        print(f"Sum of y components: {y_sum}")
        print(f"Sum of z components: {z_sum}")

        final_x_str = f"({x_sum})"
        final_y_str = f"({y_sum})"
        final_z_str = f"({z_sum})"

        print(f"Final check: Is {final_x_str} == 0 AND {final_y_str} == 0 AND {final_z_str} == 0?")

        if x_sum == 0 and y_sum == 0 and z_sum == 0:
            print("YES")
        else:
            print("NO")

    except Exception as e:
        print(f"An error occurred: {e}")

solve()
