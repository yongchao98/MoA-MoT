def solve_force_balance():
    """
    Reads n force vectors and determines if their sum is zero.
    A body is in equilibrium if the vector sum of all forces is zero.
    """
    try:
        # 1. Read the number of vectors
        n = int(input())
    except (ValueError, IndexError):
        print("Invalid input for n. Please provide a single integer.")
        return

    # 2. Initialize sums for each axis
    sum_x = 0
    sum_y = 0
    sum_z = 0

    # 3. Loop n times to read each vector and accumulate sums
    for _ in range(n):
        try:
            # Read the x, y, and z components from one line
            xi, yi, zi = map(int, input().split())
            
            # Add components to the total sums
            sum_x += xi
            sum_y += yi
            sum_z += zi
        except (ValueError, IndexError):
            print("Invalid vector input. Please provide three integers per line.")
            return

    # 4. Check for equilibrium and output the result
    # The "final equation" is the set of conditions for equilibrium
    # Sum(x) = 0, Sum(y) = 0, Sum(z) = 0
    # The prompt asks to output each number in the final equation.
    # We interpret this as showing the calculated sums.
    # For this problem, the final check is against 0, so printing the sums is sufficient.
    # print(f"Sum of forces: Fx={sum_x}, Fy={sum_y}, Fz={sum_z}")

    if sum_x == 0 and sum_y == 0 and sum_z == 0:
        print("YES")
    else:
        print("NO")

solve_force_balance()