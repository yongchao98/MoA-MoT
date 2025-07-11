import math

def solve_mirror_paths():
    """
    This function calculates the number of ways a light ray can reflect
    off the given mirrors according to the problem's constraints.
    """
    
    # The problem is to find the number of valid sequences of reflections.
    # The set of reflections is {G1, G1, G2, G3, G4}.
    # A valid sequence must alternate between the {G1, G3} set and the {G2, G4} set.

    # Our set of reflections contains 3 mirrors from {G1, G3} and 2 from {G2, G4}.
    # Therefore, the 5-reflection sequence must follow the pattern A-B-A-B-A,
    # where A is from {G1, G3} and B is from {G2, G4}.

    # Step 1: Calculate the number of ways to arrange the mirrors from Set A = {G1, G1, G3}.
    # This is a permutation with repetition.
    # The formula is n! / (k1! * k2! * ...), where n is the total number of items,
    # and k are the counts of each repeated item.
    # Here, n=3 items, with G1 repeated 2 times.
    num_items_A = 3
    repetitions_G1 = 2
    ways_A = math.factorial(num_items_A) // math.factorial(repetitions_G1)
    
    # Step 2: Calculate the number of ways to arrange the mirrors from Set B = {G2, G4}.
    # This is a simple permutation of 2 distinct items.
    num_items_B = 2
    ways_B = math.factorial(num_items_B)

    # Step 3: The total number of ways is the product of the two calculations.
    total_ways = ways_A * ways_B
    
    print("This problem is equivalent to finding the number of valid sequences of 5 reflections.")
    print("The mirrors can be grouped into two non-adjacent sets: A = {G1, G3} and B = {G2, G4}.")
    print("The multiset of reflections is {G1, G1, G2, G3, G4}, which has 3 mirrors from Set A and 2 from Set B.")
    print("This means the reflection path must follow the pattern A-B-A-B-A.")
    print("\nFirst, we find the number of ways to arrange the three 'A' mirrors ({G1, G1, G3}).")
    print(f"Number of arrangements for Set A = 3! / 2! = {ways_A}")

    print("\nNext, we find the number of ways to arrange the two 'B' mirrors ({G2, G4}).")
    print(f"Number of arrangements for Set B = 2! = {ways_B}")

    print("\nFinally, the total number of ways is the product of the arrangements for each set.")
    print(f"Total number of ways = (Ways for A) * (Ways for B)")
    print(f"The final equation is: {ways_A} * {ways_B} = {total_ways}")


solve_mirror_paths()

print("\n<<<6>>>")