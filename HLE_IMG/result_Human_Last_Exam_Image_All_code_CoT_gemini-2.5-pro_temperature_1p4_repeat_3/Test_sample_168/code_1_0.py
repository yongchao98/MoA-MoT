import math

def solve_reflection_paths():
    """
    Calculates the number of possible light ray paths based on combinatorial principles.
    """
    # Step 1: Define the number of reflections on each mirror.
    # We are given m=2 reflections on G1, and 1 on all others.
    num_G1 = 2
    num_G2 = 1
    num_G3 = 1
    num_G4 = 1

    # Step 2: Classify mirrors into Vertical (V) and Horizontal (H) groups
    # and count the number of mirrors in each group for our reflection set.
    # Vertical mirrors: G1, G3. Horizontal mirrors: G2, G4.
    vertical_mirrors_count = num_G1 + num_G3
    horizontal_mirrors_count = num_G2 + num_G4

    print("The problem is to find the number of valid sequences for 5 reflections.")
    print("The set of mirrors to be used is {G1, G1, G2, G3, G4}.")
    print("Physical constraints require the reflection sequence to alternate between Vertical (G1, G3) and Horizontal (G2, G4) mirrors.")
    print(f"\nOur set has {vertical_mirrors_count} vertical mirrors ({num_G1} G1's, {num_G3} G3) and {horizontal_mirrors_count} horizontal mirrors ({num_G2} G2, {num_G4} G4).")
    print("For a path of 5 reflections, the structure must be V-H-V-H-V to use 3 vertical and 2 horizontal mirrors.")

    # Step 3: Calculate the number of ways to arrange the vertical mirrors.
    # We are arranging the multiset {G1, G1, G3}.
    # The number of permutations is n! / (count(G1)! * count(G3)!).
    ways_V = math.factorial(vertical_mirrors_count) / math.factorial(num_G1)
    ways_V = int(ways_V)

    print(f"\nFirst, we calculate the number of ways to arrange the {vertical_mirrors_count} vertical mirrors in the three V-slots.")
    print(f"The calculation for permutations with repetition is: {vertical_mirrors_count}! / {num_G1}! = {math.factorial(vertical_mirrors_count)} / {math.factorial(num_G1)} = {ways_V}")


    # Step 4: Calculate the number of ways to arrange the horizontal mirrors.
    # We are arranging the set {G2, G4}.
    # The number of permutations is n!.
    ways_H = math.factorial(horizontal_mirrors_count)
    ways_H = int(ways_H)
    
    print(f"\nNext, we calculate the number of ways to arrange the {horizontal_mirrors_count} horizontal mirrors in the two H-slots.")
    print(f"The calculation for permutations is: {horizontal_mirrors_count}! = {math.factorial(horizontal_mirrors_count)} = {ways_H}")

    # Step 5: Calculate the total number of unique paths.
    total_ways = ways_V * ways_H
    
    print("\nThe total number of possible paths is the product of the arrangements for vertical and horizontal mirrors.")
    print(f"Final Equation: Total Ways = (Ways for Vertical Mirrors) * (Ways for Horizontal Mirrors)")
    print(f"Final Calculation: {ways_V} * {ways_H} = {total_ways}")

solve_reflection_paths()
<<<6>>>