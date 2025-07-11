import math
from collections import Counter

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def lcm(a, b):
    """Computes the least common multiple of a and b."""
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // gcd(a, b)

def lcm_list(numbers):
    """Computes the least common multiple of a list of numbers."""
    if not numbers:
        return 1
    result = numbers[0]
    for i in range(1, len(numbers)):
        result = lcm(result, numbers[i])
    return result

def get_partitions(n):
    """Generates all integer partitions of n."""
    a = [0] * (n + 1)
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[:k + 1]

def solve():
    """
    Calculates the cardinality of the coclassifier set for n=9.
    The cardinality is given by the sum over all permutations g in S_n of n!/ord(g).
    This is calculated by summing over all conjugacy classes (partitions of n).
    Formula: sum_{lambda |- n} |C_lambda| * (n! / lcm(lambda))
    """
    n = 9
    n_factorial = math.factorial(n)
    total_cardinality = 0

    print(f"Calculating the cardinality for n={n}")
    print(f"n! = {n_factorial}")
    print("-" * 30)

    partitions = get_partitions(n)
    
    for p in partitions:
        # For a partition, calculate the size of the centralizer |Z|
        counts = Counter(p)
        centralizer_size = 1
        for i, c_i in counts.items():
            centralizer_size *= math.factorial(c_i) * (i ** c_i)
        
        # Calculate the size of the conjugacy class |C|
        class_size = n_factorial // centralizer_size
        
        # Calculate the order of elements in this class (lcm of parts)
        order = lcm_list(p)
        
        # Add the summand to the total
        summand = class_size * (n_factorial // order)
        total_cardinality += summand
        
        print(f"Partition: {str(p):<20} |C|: {class_size:<7} Order: {order:<2} Summand: {class_size} * ({n_factorial} // {order}) = {summand}")

    print("-" * 30)
    print(f"Total cardinality: {total_cardinality}")

solve()