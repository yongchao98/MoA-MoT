import math
from collections import Counter

def count_permutations(numbers):
    """
    Calculates the number of distinct permutations of a list of numbers using the formula:
    n! / (n1! * n2! * ... * nk!), where n is the total number of items, and n1, n2, ...
    are the counts of each unique item.
    """
    n = len(numbers)
    counts = Counter(numbers)
    denominator = 1
    for count in counts.values():
        denominator *= math.factorial(count)
    return math.factorial(n) // denominator

def find_solutions():
    """
    This function finds all non-negative integer solutions to the Diophantine equation
    x1^2 + x2^2 + x3^2 + x4^2 + x5^2 = 2024.
    It iterates through possible values for x1, x2, x3, x4 with an ordering
    constraint (x1 <= x2 <= x3 <= x4 <= x5) to find unique sets of solutions efficiently.
    For each unique set, it calculates the number of ordered solutions (permutations)
    and adds it to the total count.
    """
    target = 2024
    limit = int(math.sqrt(target))
    
    # Pre-calculating squares and storing them in a set allows for O(1) average time complexity
    # checks to see if a number is a perfect square.
    squares_set = {i*i for i in range(limit + 1)}
    
    total_solutions = 0
    
    print(f"Finding non-negative integer solutions for x1^2 + x2^2 + x3^2 + x4^2 + x5^2 = {target}\n")
    print("Base solution sets and their contributions:")
    
    # We enforce x1 <= x2 <= x3 <= x4 <= x5 to find each unique set of numbers only once.
    # The loop limits are optimized based on this ordering, which significantly reduces the search space.
    
    # Bound for x1: 5 * x1^2 <= target  => x1 <= sqrt(target/5)
    limit1 = int(math.sqrt(target / 5))
    for x1 in range(limit1 + 1):
        s1 = x1 * x1
        
        # Bound for x2: s1 + 4 * x2^2 <= target => x2 <= sqrt((target - s1) / 4)
        limit2 = int(math.sqrt((target - s1) / 4))
        for x2 in range(x1, limit2 + 1):
            s2 = x2 * x2
            
            # Bound for x3: s1 + s2 + 3 * x3^2 <= target => x3 <= sqrt((target - s1 - s2) / 3)
            limit3 = int(math.sqrt((target - s1 - s2) / 3))
            for x3 in range(x2, limit3 + 1):
                s3 = x3 * x3
                
                # Bound for x4: s1 + s2 + s3 + 2 * x4^2 <= target => x4 <= sqrt((target - s1 - s2 - s3) / 2)
                limit4 = int(math.sqrt((target - s1 - s2 - s3) / 2))
                for x4 in range(x3, limit4 + 1):
                    s4 = x4 * x4
                    
                    current_sum = s1 + s2 + s3 + s4
                    remaining_square = target - current_sum
                    
                    if remaining_square in squares_set:
                        x5 = int(math.sqrt(remaining_square))
                        
                        # The last check to ensure the order is maintained.
                        if x5 >= x4:
                            solution_tuple = (x1, x2, x3, x4, x5)
                            perms = count_permutations(solution_tuple)
                            total_solutions += perms
                            
                            # Output each number in the equation for each unique set found.
                            num_strs = [f"{x}^2" for x in solution_tuple]
                            val_strs = [f"{x*x}" for x in solution_tuple]
                            print(f"{' + '.join(num_strs)} = {' + '.join(val_strs)} = {target}. This set gives {perms} solutions.")

    print(f"\nThe total number of non-negative integer solutions is {total_solutions}.")

if __name__ == '__main__':
    find_solutions()
    print("<<<904>>>")