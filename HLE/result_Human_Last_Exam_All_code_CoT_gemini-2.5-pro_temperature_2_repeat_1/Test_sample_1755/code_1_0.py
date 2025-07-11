import sys

def solve():
    """
    Reads force vectors from stdin, calculates their sum,
    and determines if they are in equilibrium.
    """
    try:
        # The first line of input is the number of vectors, n.
        n_line = sys.stdin.readline()
        if not n_line.strip():
            # Handle empty input gracefully.
            print("NO")
            return
            
        n = int(n_line)
        if not (1 <= n <= 100):
            # As per problem constraints, n is between 1 and 100.
            # Handle invalid n, although not specified in the problem output.
            print("NO")
            return
            
        vectors = []
        for _ in range(n):
            line = sys.stdin.readline()
            vectors.append(list(map(int, line.split())))

        # Separate coordinates for easier summation and display
        x_coords = [v[0] for v in vectors]
        y_coords = [v[1] for v in vectors]
        z_coords = [v[2] for v in vectors]

        sum_x = sum(x_coords)
        sum_y = sum(y_coords)
        sum_z = sum(z_coords)

        # Output each number in the final equation as requested.
        # The .replace() is for cleaner display of negative numbers.
        print(' + '.join(map(str, x_coords)).replace('+ -', '- ') + f' = {sum_x}')
        print(' + '.join(map(str, y_coords)).replace('+ -', '- ') + f' = {sum_y}')
        print(' + '.join(map(str, z_coords)).replace('+ -', '- ') + f' = {sum_z}')

        # Check for equilibrium and print the final result.
        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            print("NO")

    except (ValueError, IndexError, IOError):
        # If input is malformed or missing, the body is not in equilibrium.
        print("NO")

solve()