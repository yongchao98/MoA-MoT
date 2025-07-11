import sys

def solve():
    """
    Determines if a body is in equilibrium based on a set of force vectors.
    """
    try:
        # 1. Read the number of vectors.
        n = int(sys.stdin.readline())

        # 2. Initialize sum of forces for x, y, z axes.
        force_sum = [0, 0, 0]

        # 3. Loop n times to read all vectors.
        for _ in range(n):
            # 4. Read the three components of a force vector.
            line = sys.stdin.readline()
            if not line:
                continue
            x, y, z = map(int, line.strip().split())
            
            # 5. Add components to the total sum.
            force_sum[0] += x
            force_sum[1] += y
            force_sum[2] += z

        # 7. Print the final resultant force vector.
        print(f"Resultant Force Vector = ({force_sum[0]}, {force_sum[1]}, {force_sum[2]})")
        
        # 8. Check for equilibrium and print the result.
        if force_sum[0] == 0 and force_sum[1] == 0 and force_sum[2] == 0:
            print("YES")
        else:
            print("NO")

    except (IOError, ValueError) as e:
        print(f"An error occurred: {e}", file=sys.stderr)

solve()