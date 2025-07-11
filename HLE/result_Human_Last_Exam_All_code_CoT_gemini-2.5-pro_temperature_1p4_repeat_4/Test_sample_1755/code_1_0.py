import sys

def solve():
    """
    Reads n vectors, calculates their sum, and determines if a body is in equilibrium.
    """
    try:
        # Read the number of vectors
        n_str = sys.stdin.readline()
        if not n_str:
            print("NO")
            return
        n = int(n_str)
        
        if not (1 <= n <= 100):
            # Handle edge case of invalid n, though problem statement guarantees 1 <= n <= 100
            print("NO")
            return

        all_coords = []
        for _ in range(n):
            line = sys.stdin.readline()
            if not line:
                break
            all_coords.append(list(map(int, line.split())))

        if len(all_coords) != n:
            # Handle case of incomplete input
            print("NO")
            return
            
        x_coords = [vec[0] for vec in all_coords]
        y_coords = [vec[1] for vec in all_coords]
        z_coords = [vec[2] for vec in all_coords]

        sum_x = sum(x_coords)
        sum_y = sum(y_coords)
        sum_z = sum(z_coords)

        # Format the summation equation strings
        x_eq_str = " + ".join(f"({i})" if i < 0 else str(i) for i in x_coords)
        y_eq_str = " + ".join(f"({i})" if i < 0 else str(i) for i in y_coords)
        z_eq_str = " + ".join(f"({i})" if i < 0 else str(i) for i in z_coords)

        # Print the final equations
        print(f"x: {x_eq_str} = {sum_x}")
        print(f"y: {y_eq_str} = {sum_y}")
        print(f"z: {z_eq_str} = {sum_z}")
        
        # Determine and print the final result
        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            print("NO")

    except (IOError, ValueError) as e:
        # Handle potential errors during input reading/conversion
        print(f"An error occurred: {e}")
        print("NO")

solve()