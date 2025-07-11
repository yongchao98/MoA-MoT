import math

def solve_tiling_problem():
    """
    Finds the area of the smallest integer-length rectangle admitting a
    non-guillotine tiling by squares from a given set.
    """
    # Step 1: Define the set of allowed square side lengths.
    S = [2, 3, 5, 7]
    print(f"The given set of allowed square side lengths is S = {S}.")
    print("-" * 50)

    # Step 2: Explain the condition for a non-guillotine tiling.
    print("A tiling 'not constructable with glass-cuts' (or non-guillotine) is one that cannot be divided")
    print("into two pieces by a single straight cut running from one edge of the rectangle to the opposite edge.")
    print("\nThe smallest and most basic structure that creates such a tiling with squares is a 'pinwheel' pattern.")
    print("This pattern involves four squares meeting at a single interior point.")
    print("-" * 50)

    # Step 3: Describe the construction of the minimal pinwheel.
    print("For four squares to form a pinwheel shape that tiles a rectangle, they must consist of two pairs")
    print("of identical squares. Let's say we have two squares of side 'a' and two of side 'b'.")
    print("When arranged in a pinwheel, they perfectly tile a larger square with a side length of (a + b).")
    print("-" * 50)

    # Step 4: Find the smallest possible pinwheel dimensions from the set S.
    print(f"To find the smallest rectangle with a non-guillotine tiling, we must find the smallest possible")
    print(f"pinwheel square. This means finding the minimum side length of '(a + b)' where 'a' and 'b'")
    print(f"are distinct side lengths from the set S = {S}.")

    min_side_sum = float('inf')
    best_a = -1
    best_b = -1

    # Step 5: Iterate through all unique pairs from S to find the minimum sum.
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            a = S[i]
            b = S[j]
            current_sum = a + b
            if current_sum < min_side_sum:
                min_side_sum = current_sum
                best_a = a
                best_b = b

    print(f"\nThe smallest sum of two distinct side lengths from S is found with sides {best_a} and {best_b}.")
    
    # Step 6: Identify the smallest rectangle and its properties.
    smallest_side = best_a + best_b
    print(f"This results in a pinwheel that forms a square of side {best_a} + {best_b} = {smallest_side}.")
    print(f"\nTherefore, the smallest rectangle that meets the criteria is a {smallest_side}x{smallest_side} square.")
    print("This rectangle has a non-guillotine tiling (two {0}x{0} squares and two {1}x{1} squares).".format(best_a, best_b))
    print("It also can be tiled in a guillotine fashion (with a single {0}x{0} square), as {0} is in the set S.".format(smallest_side))
    print("-" * 50)

    # Step 7: Calculate and display the final area.
    length = smallest_side
    width = smallest_side
    area = length * width

    print("The area of this smallest rectangle is calculated as follows:")
    print(f"{length} * {width} = {area}")


# Run the solver
solve_tiling_problem()
<<<25>>>