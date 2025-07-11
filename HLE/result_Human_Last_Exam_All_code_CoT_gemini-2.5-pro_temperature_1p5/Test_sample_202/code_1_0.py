import sys
import math

def solve_triangle_distance():
    """
    Calculates the shortest distance between the maximum and minimum value nodes
    in a triangulated grid based on the values at vertices A, B, and C.
    """
    # --- Read input values for a, b, c ---
    # The script expects three numeric arguments from the command line for a, b, and c.
    if len(sys.argv) != 4:
        print("Usage: python your_script_name.py <a> <b> <c>")
        print("Please provide the numeric values for the vertices A, B, and C.")
        return

    try:
        a = float(sys.argv[1])
        b = float(sys.argv[2])
        c = float(sys.argv[3])
    except ValueError:
        print("Error: The values for a, b, and c must be numeric.")
        return

    # --- Determine the case based on the equality of a, b, c ---
    
    # Case 1: a, b, and c are all equal
    if a == b and b == c:
        distance = 0
        print(f"All vertices have the same value ({a}).")
        print("Therefore, all nodes in the triangle have this value.")
        print("The shortest distance between a point with the maximum value and a point with the minimum value is 0.")

    # Case 2: Exactly two of a, b, c are equal
    elif a == b or b == c or a == c:
        distance = math.sqrt(3) / 2
        vals = {'A': a, 'B': b, 'C': c}
        val_list = sorted(vals.values())
        min_val = val_list[0]
        max_val = val_list[2]

        if val_list[0] == val_list[1]: # One max, two mins
            max_vertex = [k for k, v in vals.items() if v == max_val][0]
            min_side_vertices = [k for k, v in vals.items() if v == min_val]
            side_name = "".join(sorted(min_side_vertices))
            print(f"Maximum value ({max_val}) is at vertex {max_vertex}.")
            print(f"Minimum value ({min_val}) occurs along the side {side_name}.")
        else: # Two maxes, one min
            min_vertex = [k for k, v in vals.items() if v == min_val][0]
            max_side_vertices = [k for k, v in vals.items() if v == max_val]
            side_name = "".join(sorted(max_side_vertices))
            print(f"Maximum value ({max_val}) occurs along the side {side_name}.")
            print(f"Minimum value ({min_val}) is at vertex {min_vertex}.")

        print(f"The shortest distance is the altitude of the triangle, sqrt(3)/2, which is approximately {distance:.5f}.")

    # Case 3: a, b, c are all distinct
    else:
        distance = 1
        vals = {'A': a, 'B': b, 'C': c}
        min_val = min(vals.values())
        max_val = max(vals.values())
        min_vertex = [k for k, v in vals.items() if v == min_val][0]
        max_vertex = [k for k, v in vals.items() if v == max_val][0]
        
        print(f"The values a, b, c are all distinct.")
        print(f"Maximum value ({max_val}) is at vertex {max_vertex}.")
        print(f"Minimum value ({min_val}) is at vertex {min_vertex}.")
        print(f"The shortest distance is the side length of the triangle, which is 1.")

if __name__ == '__main__':
    solve_triangle_distance()