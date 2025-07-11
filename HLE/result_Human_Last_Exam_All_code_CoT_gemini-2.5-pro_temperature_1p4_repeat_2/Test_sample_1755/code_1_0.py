import sys

def solve():
    """
    Reads n vectors, calculates the sum of forces for each coordinate,
    and determines if the body is in equilibrium.
    """
    try:
        n_str = sys.stdin.readline()
        if not n_str:
            print("Error: Input is empty.")
            return
        n = int(n_str)
        if not (1 <= n <= 100):
            # This check is based on problem constraints, optional here
            # print("Error: n is out of range (1-100).")
            pass

        xs, ys, zs = [], [], []
        for _ in range(n):
            line = sys.stdin.readline()
            if not line:
                break
            parts = list(map(int, line.split()))
            if len(parts) == 3:
                xs.append(parts[0])
                ys.append(parts[1])
                zs.append(parts[2])

        sum_x = sum(xs)
        sum_y = sum(ys)
        sum_z = sum(zs)

        # Build and print the equation strings
        x_eq = " + ".join(map(str, xs)).replace("+ -", "- ")
        y_eq = " + ".join(map(str, ys)).replace("+ -", "- ")
        z_eq = " + ".join(map(str, zs)).replace("+ -", "- ")
        
        print(f"x: {x_eq} = {sum_x}")
        print(f"y: {y_eq} = {sum_y}")
        print(f"z: {z_eq} = {sum_z}")

        if sum_x == 0 and sum_y == 0 and sum_z == 0:
            print("YES")
        else:
            print("NO")

    except (ValueError, IndexError):
        print("Error: Invalid input format.")

solve()