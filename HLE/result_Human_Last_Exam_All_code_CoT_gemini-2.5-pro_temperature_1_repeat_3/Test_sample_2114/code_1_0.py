import heapq
import math

def frobenius_number(numbers):
    """
    Calculates the Frobenius number for a set of integers with gcd=1.
    The algorithm is based on finding the shortest path in a residue graph
    using Dijkstra's algorithm.
    """
    # Remove duplicates and numbers that are combinations of others
    # A simple sort and check for divisibility helps, but isn't a full reduction
    numbers = sorted(list(set(numbers)))
    
    # Check for GCD
    if len(numbers) > 1:
        current_gcd = numbers[0]
        for i in range(1, len(numbers)):
            current_gcd = math.gcd(current_gcd, numbers[i])
        if current_gcd > 1:
            # The Frobenius number is not defined (or infinite) if gcd > 1.
            # The problem implies a valid set.
            return -1 # Or raise an error
    elif len(numbers) == 1 and numbers[0] > 1:
        return -1


    # The problem reduces to {2,3} since 9 is a multiple of 3.
    # For two variables {a,b}, the formula is ab-a-b.
    if len(numbers) == 2:
        return numbers[0] * numbers[1] - numbers[0] - numbers[1]

    # General algorithm for n >= 2
    a1 = numbers[0]
    other_numbers = numbers[1:]
    
    # dists[i] will store the smallest number congruent to i (mod a1)
    # that can be formed by a sum from the set.
    dists = {i: float('inf') for i in range(a1)}
    dists[0] = 0
    
    # Priority queue for Dijkstra's algorithm
    pq = [(0, 0)] # (distance, vertex)
    
    while pq:
        d, u = heapq.heappop(pq)
        
        if d > dists[u]:
            continue
            
        for num in other_numbers:
            v = (u + num) % a1
            new_dist = d + num
            if new_dist < dists[v]:
                dists[v] = new_dist
                heapq.heappush(pq, (new_dist, v))
                
    max_dist = 0
    for dist in dists.values():
        if dist == float('inf'):
            # Should not happen if gcd is 1
            return -1
        if dist > max_dist:
            max_dist = dist
            
    return max_dist - a1

# Based on the plan, the values are derived as:
X1 = 4
X2 = 2
X3 = 3

# Calculate the set for the Frobenius number problem
a1 = math.ceil(X1 + X2 + X3)
a2 = math.ceil(X2)
a3 = math.ceil(X3)

number_set = [a1, a2, a3]

# Calculate the Frobenius number
g = frobenius_number(number_set)

# As requested, output the numbers in the final equation.
# Here we show the set and its calculated Frobenius number.
# The set is sorted for clarity.
sorted_set = sorted(number_set)
print(f"The set of integers is {{{sorted_set[0]}, {sorted_set[1]}, {sorted_set[2]}}}.")
print(f"The Frobenius Number is g({sorted_set[0]}, {sorted_set[1]}, {sorted_set[2]}) = {g}")
