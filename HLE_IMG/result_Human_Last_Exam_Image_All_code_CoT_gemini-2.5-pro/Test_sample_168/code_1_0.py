import math

def solve_light_path_ways():
    """
    Calculates the number of ways a light ray can travel given a specific
    number of reflections on a set of four mirrors.
    """
    # The problem specifies the number of reflections on each mirror.
    # Total reflections (n) is the sum of reflections on each mirror.
    # m=2 reflections on G1.
    # 1 reflection on G2, 1 on G3, 1 on G4.
    
    # These are the counts of identical items in our multiset of reflections.
    reflections_G1 = 2
    reflections_G2 = 1
    reflections_G3 = 1
    reflections_G4 = 1

    # n is the total number of reflections (the size of the multiset).
    n = reflections_G1 + reflections_G2 + reflections_G3 + reflections_G4

    # Calculate the factorials needed for the permutation formula.
    n_factorial = math.factorial(n)
    k1_factorial = math.factorial(reflections_G1)
    k2_factorial = math.factorial(reflections_G2)
    k3_factorial = math.factorial(reflections_G3)
    k4_factorial = math.factorial(reflections_G4)

    # The denominator is the product of the factorials of the counts of each identical item.
    denominator = k1_factorial * k2_factorial * k3_factorial * k4_factorial

    # The number of ways is the result of the permutation formula.
    number_of_ways = n_factorial // denominator

    # Print the step-by-step explanation and calculation.
    print("The problem asks for the number of ways a light ray can be drawn, which corresponds to the number of unique sequences of reflections.")
    print("This can be solved using the formula for permutations of a multiset.")
    print("")
    print(f"Total number of reflections (n) = {reflections_G1} + {reflections_G2} + {reflections_G3} + {reflections_G4} = {n}")
    print(f"The multiset of reflections is {{G1, G1, G2, G3, G4}}.")
    print("The counts for each type of reflection are:")
    print(f"  - Number of reflections on G1 (k1) = {reflections_G1}")
    print(f"  - Number of reflections on G2 (k2) = {reflections_G2}")
    print(f"  - Number of reflections on G3 (k3) = {reflections_G3}")
    print(f"  - Number of reflections on G4 (k4) = {reflections_G4}")
    print("")
    print("The formula is: n! / (k1! * k2! * k3! * k4!)")
    print("Substituting the values into the equation:")
    print(f"  {n}! / ({reflections_G1}! * {reflections_G2}! * {reflections_G3}! * {reflections_G4}!)")
    print(f"= {n_factorial} / ({k1_factorial} * {k2_factorial} * {k3_factorial} * {k4_factorial})")
    print(f"= {n_factorial} / {denominator}")
    print(f"= {number_of_ways}")
    print("")
    print(f"So, there are {number_of_ways} ways to draw the light ray path.")

solve_light_path_ways()
<<<60>>>