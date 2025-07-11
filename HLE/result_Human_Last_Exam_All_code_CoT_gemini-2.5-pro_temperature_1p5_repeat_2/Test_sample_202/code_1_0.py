import math

def solve():
    """
    Reads input values for n, a, b, c, calculates the shortest distance
    between the max and min valued nodes, and prints the result.
    """
    try:
        # Prompt the user for input and read it from a single line.
        # Example input: 3 10 10 4
        print("Please enter n, a, b, c separated by spaces:")
        input_line = input()
        parts = input_line.split()

        if len(parts) < 4:
            raise ValueError("Four values are required.")

        n = int(parts[0])
        a = float(parts[1])
        b = float(parts[2])
        c = float(parts[3])

        if n <= 0:
            raise ValueError("n must be a positive integer.")

        vals = [a, b, c]
        max_val = max(vals)
        min_val = min(vals)

        # Case 1: All vertex values are the same.
        # All nodes have the same value, so the distance is 0.
        if max_val == min_val:
            print("The values at all vertices are the same.")
            print("Distance = 0.0")
            return # Added return to avoid printing the answer tag on a separate line

        max_count = vals.count(max_val)
        min_count = vals.count(min_val)
        
        # Case 2: Max and Min values are at unique vertices.
        # This occurs when a > b > c (or any permutation).
        if max_count == 1 and min_count == 1:
            print("The maximum and minimum values are at unique vertices.")
            print("Distance = 1.0")

        # Case 3: One extremum is on an edge, the other is at a vertex.
        # This occurs when a=b > c or a > b=c (or any permutation).
        else:
            if n % 2 == 0:
                distance = math.sqrt(3) / 2
                print(f"One extremum is on a side, the other at a vertex, and n ({n}) is even.")
                print(f"Distance = sqrt(3)/2 = {distance:.12f}")
            else:
                distance = math.sqrt((1 / (2 * n))**2 + 3 / 4)
                print(f"One extremum is on a side, the other at a vertex, and n ({n}) is odd.")
                # The ^ symbol is for display; Python uses ** for exponentiation.
                print(f"Distance = sqrt((1 / (2 * {n}))^2 + 3/4) = {distance:.12f}")

    except (ValueError, IndexError) as e:
        print(f"Invalid input: {e}. Please enter four numbers separated by spaces (e.g., '3 10 10 4').")


# Execute the solution
solve()