import sys

def check_equilibrium():
    """
    Reads n force vectors and determines if their sum is zero.
    A body is in equilibrium if the net force vector is (0, 0, 0).
    """
    try:
        # 1. Read the number of vectors.
        n = int(sys.stdin.readline())
    except (ValueError, IndexError):
        # In case of empty or invalid input for n.
        print("NO")
        return

    # 2. Initialize sums of vector components.
    sum_x, sum_y, sum_z = 0, 0, 0

    # 3. Loop n times to read each vector.
    for _ in range(n):
        try:
            # 4. Read the three components of a force vector.
            # The map(int, ...) function efficiently converts them to integers.
            xi, yi, zi = map(int, sys.stdin.readline().split())

            # 5. Add components to their respective sums.
            sum_x += xi
            sum_y += yi
            sum_z += zi
        except (ValueError, IndexError):
            # This handles cases of malformed input lines. For this problem,
            # we assume valid input, but this makes the code more robust.
            # If an invalid line is found, the body cannot be guaranteed
            # to be in equilibrium, so we can default to "NO".
            print("NO")
            return

    # 6. Check if the net force is the zero vector.
    if sum_x == 0 and sum_y == 0 and sum_z == 0:
        # 7. Print "YES" for equilibrium.
        print("YES")
    else:
        # 7. Print "NO" otherwise.
        print("NO")

check_equilibrium()