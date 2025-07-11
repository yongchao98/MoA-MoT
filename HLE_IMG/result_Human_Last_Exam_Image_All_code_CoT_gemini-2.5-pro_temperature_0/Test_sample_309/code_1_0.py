import math

def gcd(a, b):
    """Computes the greatest common divisor of two integers."""
    while b:
        a, b = b, a % b
    return a

def frobenius_number(numbers):
    """
    Computes the Frobenius number for a set of integers.
    This implementation uses a common algorithm for three numbers.
    It can be generalized but is shown here for the specific case.
    """
    # Check if the gcd of the set is 1
    if len(numbers) < 2:
        raise ValueError("Need at least two numbers.")
    
    result_gcd = numbers[0]
    for i in range(1, len(numbers)):
        result_gcd = gcd(result_gcd, numbers[i])
    
    if result_gcd != 1:
        # The Frobenius number is not defined (infinite) if gcd is not 1.
        # However, the problem implies a solution exists.
        # This check is for correctness.
        print(f"The GCD of the set {numbers} is {result_gcd}, not 1. The Frobenius number is technically infinite.")
        # For this problem, we proceed as the context implies coprimality for the overall problem.
        pass

    # Sort numbers for easier processing
    numbers.sort()
    a1 = numbers[0]
    
    # Using the algorithm based on congruences and shortest paths in a cyclic graph
    # For each residue r mod a1, find the smallest number t_r representable by the rest of the numbers
    # t_r = min{a2*y + a3*z + ... | a2*y + a3*z + ... = r (mod a1)}
    
    # For our set {9, 51, 53}, a1 = 9. We need to find t_r for r in {1..8}
    # using numbers {51, 53}.
    
    # dist[r] will store the smallest number = r (mod a1)
    dist = {i: float('inf') for i in range(a1)}
    dist[0] = 0
    
    # Using Dijkstra-like approach
    import heapq
    pq = [(0, 0)] # (distance, residue)
    
    while pq:
        d, u = heapq.heappop(pq)
        
        if d > dist[u]:
            continue
        
        for i in range(1, len(numbers)):
            num = numbers[i]
            v = (u + num) % a1
            new_dist = d + num
            
            if new_dist < dist[v]:
                dist[v] = new_dist
                heapq.heappush(pq, (new_dist, v))
                
    # The Frobenius number is the maximum of these smallest numbers, minus a1
    max_t_r = 0
    for r in range(a1):
        if dist[r] == float('inf'):
            # This case shouldn't happen if gcd is 1
            return -1 # Indicates an error
        if dist[r] > max_t_r:
            max_t_r = dist[r]
            
    frobenius = max_t_r - a1
    return frobenius

# The set of numbers determined from the problem analysis
number_set = [9, 51, 53]

# Calculate the Frobenius number
frobenius_result = frobenius_number(number_set)

# The problem asks to show the numbers in the final equation.
# The Frobenius number g is the result of g({a, b, c}).
a, b, c = sorted(number_set)
print(f"The set of numbers for the Frobenius calculation is {{{a}, {b}, {c}}}.")
print(f"The Frobenius number g({a}, {b}, {c}) is the largest integer that cannot be expressed in the form {a}*x + {b}*y + {c}*z for non-negative integers x, y, z.")
print(f"The calculated Frobenius number is: {frobenius_result}")
