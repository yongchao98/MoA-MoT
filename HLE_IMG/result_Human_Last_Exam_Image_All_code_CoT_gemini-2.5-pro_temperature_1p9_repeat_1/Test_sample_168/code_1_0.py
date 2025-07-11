import math

def solve_reflection_paths():
    """
    Calculates the number of possible light ray paths based on the number of reflections.
    """
    # The problem specifies the number of reflections on each mirror.
    # Number of reflections on mirror G1.
    m = 2
    # Number of reflections on mirrors G2, G3, and G4 is one each.
    n_g2 = 1
    n_g3 = 1
    n_g4 = 1

    # The total number of reflections is the sum of reflections on all mirrors.
    # This determines the length of the reflection sequence.
    total_reflections = m + n_g2 + n_g3 + n_g4

    # The problem asks for the number of unique ways to arrange the sequence of reflections.
    # The set of reflections is {G1, G1, G2, G3, G4}.
    # This is a permutation with repetition problem.
    # The formula is n! / (k1! * k2! * ...), where n is the total number of items,
    # and k1, k2, ... are the counts of each distinct item.
    # In our case, n = 5, k1 (for G1) = 2, k2 (for G2) = 1, k3 (for G3) = 1, k4 (for G4) = 1.
    # The formula simplifies to 5! / (2! * 1! * 1! * 1!), which is 5! / 2!.

    # Calculate the factorials needed for the formula.
    total_reflections_factorial = math.factorial(total_reflections)
    g1_reflections_factorial = math.factorial(m)
    
    # Calculate the total number of ways.
    number_of_ways = total_reflections_factorial // g1_reflections_factorial
    
    # Print the explanation and the step-by-step calculation.
    print("To find the number of ways, we calculate the number of unique sequences of reflections.")
    print(f"Number of reflections on G1 = {m}")
    print(f"Number of reflections on G2 = {n_g2}")
    print(f"Number of reflections on G3 = {n_g3}")
    print(f"Number of reflections on G4 = {n_g4}")
    print(f"Total number of reflections in the sequence = {m} + {n_g2} + {n_g3} + {n_g4} = {total_reflections}")
    print("\nThis is a permutation with repetition problem.")
    print("The number of ways is calculated using the formula: (Total reflections)! / (Repetitions of G1)!")
    
    print("\nFinal Equation:")
    print(f"{total_reflections}! / {m}! = {total_reflections_factorial} / {g1_reflections_factorial} = {number_of_ways}")
    
    # Output the final answer in the required format.
    print(f"\n<<<{number_of_ways}>>>")

solve_reflection_paths()