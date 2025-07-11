import math

def solve():
    """
    Finds the number of ordered triplets (x, y, z) of positive integers
    such that x*y*z = 216 and x^(log_2 y) = y^(log_2 z).
    
    The second condition simplifies to (y = 1) or (x = z).
    """
    
    volume = 216
    
    # Case 1: y = 1
    # We need to find the number of pairs (x, z) such that x*z = 216.
    # This is equal to the number of divisors of 216.
    count_y_is_1 = 0
    for i in range(1, volume + 1):
        if volume % i == 0:
            count_y_is_1 += 1
            
    # Case 2: x = z
    # We need to find the number of integers x such that x^2 divides 216.
    # We only need to check x up to floor(sqrt(216)).
    count_x_equals_z = 0
    limit = int(math.sqrt(volume))
    for x in range(1, limit + 1):
        if volume % (x * x) == 0:
            # y = volume / (x*x) is guaranteed to be a positive integer
            count_x_equals_z += 1
            
    # The two cases are disjoint because x^2 = 216 has no integer solution for x.
    # Total count is the sum of the counts from both cases.
    total_permutations = count_y_is_1 + count_x_equals_z
    
    # The prompt asks to output each number in the final equation.
    print(f"{count_y_is_1} + {count_x_equals_z} = {total_permutations}")

solve()
<<<20>>>