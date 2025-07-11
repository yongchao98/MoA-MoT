def solve_force_equilibrium():
    """
    Determines if a body is in equilibrium based on a set of force vectors.
    """
    try:
        # Read the number of force vectors.
        n = int(input())
        
        # Initialize sums for each coordinate to zero.
        sum_x = 0
        sum_y = 0
        sum_z = 0
        
        # Loop n times to read each force vector.
        for _ in range(n):
            # Read the line containing the three coordinates and split them.
            line = input().split()
            # Convert coordinates to integers.
            xi = int(line[0])
            yi = int(line[1])
            zi = int(line[2])
            
            # Add the coordinates to their respective sums.
            sum_x += xi
            sum_y += yi
            sum_z += zi
            
        # Check if the body is in equilibrium.
        # The condition is that the sum of forces in all three dimensions is zero.
        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            print("NO")
            
    except (ValueError, IndexError):
        # Handle potential errors with input format.
        print("Invalid input format.")

solve_force_equilibrium()