import sys

def solve():
    """
    Determines if a body is in equilibrium based on a set of force vectors.
    """
    try:
        # 1. Read the number of vectors.
        n = int(sys.stdin.readline())
        
        # 2. Initialize sum of forces to zero.
        # force_sum[0] for x, force_sum[1] for y, force_sum[2] for z
        force_sum = [0, 0, 0]
        
        # 3. Loop n times to read all vectors.
        for _ in range(n):
            # 4. Read each vector and add it to the sum.
            line = sys.stdin.readline()
            components = list(map(int, line.split()))
            force_sum[0] += components[0]
            force_sum[1] += components[1]
            force_sum[2] += components[2]
            
        # 5. Output the resultant force vector components (the final equation).
        print(f"Sum of forces: Fx={force_sum[0]}, Fy={force_sum[1]}, Fz={force_sum[2]}")
        
        # 6. Check for equilibrium and print the result.
        if force_sum[0] == 0 and force_sum[1] == 0 and force_sum[2] == 0:
            print("YES")
        else:
            print("NO")
            
    except (ValueError, IndexError):
        print("Invalid input format. Please check the input data.")

solve()