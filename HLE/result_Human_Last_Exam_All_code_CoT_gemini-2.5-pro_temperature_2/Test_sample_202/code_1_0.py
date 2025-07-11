import math

# --- User-configurable inputs ---
# n: the number of divisions on each side of the triangle
# a, b, c: the numbers placed on vertices A, B, C respectively
n = 4
a = 10
b = 20
c = 20
# --- End of inputs ---

# Store vertex values for easy analysis
vals = [a, b, c]
min_val = min(vals)
max_val = max(vals)

# Case 1: All vertex values are the same (a=b=c)
if min_val == max_val:
    distance = 0.0
    print("The values at all three vertices A, B, and C are equal.")
    print("The shortest distance is 0.")
    print("\nFinal Equation:")
    print("d = 0")
    print(f"\nResult: {distance}")

else:
    # Count how many vertices have the min and max values
    num_min = vals.count(min_val)
    num_max = vals.count(max_val)

    # Case 2: The vertex values a, b, c are all distinct
    # The max is at one vertex and the min is at another.
    if num_min == 1 and num_max == 1:
        distance = 1.0
        print("The values at vertices A, B, and C are all distinct.")
        print("The distance is between two vertices of the main triangle.")
        print("\nFinal Equation:")
        print("d = 1")
        print(f"\nResult: {distance}")

    # Case 3: Two vertices share a value (either min or max)
    # The extreme value is on a line, and the other extreme is at a point.
    else:
        # Subcase 3a: n is even
        if n % 2 == 0:
            print("Two vertices share an extreme value, and n is even.")
            print("The shortest distance is the altitude of the triangle.")
            distance = math.sqrt(3) / 2
            
            print("\nFinal Equation:")
            print("d = sqrt(3) / 2")
            print(f"d = {math.sqrt(3)} / 2")
            print(f"d = {distance}")
        
        # Subcase 3b: n is odd
        else:
            print("Two vertices share an extreme value, and n is odd.")
            print("The shortest distance is from a vertex to the closest node on the opposite side.")
            distance = math.sqrt(3 * n**2 + 1) / (2 * n)
            
            print("\nFinal Equation:")
            # Show the calculation step-by-step with the numbers
            term1 = 3 * n**2
            numerator = term1 + 1
            denominator = 2 * n
            
            print(f"d = sqrt(3 * n^2 + 1) / (2 * n)")
            print(f"d = sqrt(3 * {n}^2 + 1) / (2 * {n})")
            print(f"d = sqrt({term1} + 1) / {denominator}")
            print(f"d = sqrt({numerator}) / {denominator}")
            print(f"d = {math.sqrt(numerator)} / {denominator}")
            print(f"d = {distance}")

<<<0.8660254037844386>>>