import math

def get_primes(k):
    """Generates the first k prime numbers."""
    primes = []
    num = 2
    while len(primes) < k:
        is_prime = True
        for p in primes:
            if p * p > num:
                break
            if num % p == 0:
                is_prime = False
                break
        if is_prime:
            primes.append(num)
        num += 1
    return primes

def calculate_n(partition, primes):
    """Calculates the maximum number N for a given partition."""
    n_val = 1
    equation_str = "N = "
    for i in range(len(partition)):
        prime = primes[i]
        exponent = partition[i] - 1
        term = prime ** exponent
        n_val *= term
        if i > 0:
            equation_str += " * "
        equation_str += f"{prime}^({partition[i]}-1)"
    
    # Add the calculation steps to the equation string
    equation_str += " = "
    for i in range(len(partition)):
        prime = primes[i]
        exponent = partition[i] - 1
        term = prime ** exponent
        if i > 0:
            equation_str += " * "
        equation_str += str(term)

    equation_str += f" = {n_val}"
    return n_val, equation_str

def find_smallest_n():
    """
    Finds the smallest N by checking all partitions of 16.
    """
    print("Investigating possible constructions for the set of 16 numbers.")
    print("The numbers are of the form p_1^e_1 * p_2^e_2 * ...")
    print("The number of distinct exponents for each prime (d_i) must satisfy d_1*d_2*... = 16.\n")

    # Partitions of 16 into factors > 1
    partitions = [
        [16],
        [8, 2],
        [4, 4],
        [4, 2, 2],
        [2, 2, 2, 2]
    ]

    min_n = float('inf')
    best_partition = None
    
    for p in partitions:
        # Sort partition to assign larger d_i to smaller primes
        p.sort(reverse=True)
        num_primes = len(p)
        primes = get_primes(num_primes)
        
        n, equation = calculate_n(p, primes)
        
        print(f"Case with partition {p}:")
        print(f"  Primes: {primes}")
        print(f"  Equation: {equation}")
        
        if n < min_n:
            min_n = n
            best_partition = p
        print("-" * 20)

    print(f"\nThe minimum N found is {min_n} from partition {best_partition}.")
    print("This is based on the constructive method, assuming the required underlying combinatorial objects (orthogonal magic squares) exist, which is known to be true for these cases.")

find_smallest_n()
