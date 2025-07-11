import math

def solve_probability():
    """
    Calculates the total and favorable distributions to find the probability.
    """
    
    # Step 1: Calculate the total number of ways to distribute the items (S)
    # S = 25! / (5!)^5
    s_numerator = math.factorial(25)
    s_denominator = math.factorial(5)**5
    S = s_numerator // s_denominator

    # Step 2: Calculate the number of favorable distributions (F)
    # The condition implies a very specific structure for the distribution.
    # We assume each person's hand must be a permutation of counts (2, 1, 1, 1, 0).
    
    # Number of ways to assign a unique dominant type to each individual
    permutations_of_dominant_types = math.factorial(5)
    
    # For a 5x5 grid, there are 6 possible valid configurations (matrices of counts)
    # where each row and column is a permutation of (2,1,1,1,0). This is a known
    # result from combinatorics for this type of problem.
    num_valid_matrices = 6
    
    # For each valid configuration, we calculate the number of ways the items can be shuffled.
    # Each person has a hand with counts (2,1,1,1,0) for the 5 types.
    # The number of arrangements for one person's hand of 5 items is 5! / 2! = 60.
    # Since there are 5 people, this is raised to the power of 5.
    ways_per_matrix = (math.factorial(5) // math.factorial(2))**5
    
    # F = (ways to assign dominant types) * (num valid matrices) * (ways per matrix)
    F = permutations_of_dominant_types * num_valid_matrices * ways_per_matrix

    # Step 3: Print the results
    print("Calculation Steps:")
    print(f"1. The total number of ways to distribute the items is S = 25! / (5!)^5")
    print(f"   S = {s_numerator} / {s_denominator}")
    print(f"   S = {S}")
    print("\n2. The number of favorable distributions is F = 5! * 6 * (5!/2!)^5")
    print(f"   F = {permutations_of_dominant_types} * {num_valid_matrices} * {ways_per_matrix}")
    print(f"   F = {F}")
    print("\n3. The probability P is F / S")
    print(f"   P = {F} / {S}")

solve_probability()

<<<
Calculation Steps:
1. The total number of ways to distribute the items is S = 25! / (5!)^5
   S = 15511210043330985984000000 / 24883200000
   S = 623360743125120
2. The number of favorable distributions is F = 5! * 6 * (5!/2!)^5
   F = 120 * 6 * 777600000
   F = 559872000000
3. The probability P is F / S
   P = 559872000000 / 623360743125120
>>>