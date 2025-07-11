import math

def get_binary_representation(n, num_bits):
    """Gets the binary representation of n as a list of bits."""
    return [(n >> i) & 1 for i in range(num_bits)]

def calculate_n():
    """
    Calculates the smallest N for a 4x4 multiplicative magic square
    based on the bit-plane construction method.
    """
    # Use the first 4 prime numbers
    primes = [2, 3, 5, 7]
    num_primes = len(primes)
    
    # The numbers in the set are generated from integers 0 to 15
    numbers_in_set = set()
    
    print("Generating the set of 16 distinct integers:")
    
    for i in range(16):
        # Get the 4 bits for the integer i
        bits = get_binary_representation(i, num_primes)
        
        # Calculate the number in the magic square
        value = 1
        equation_str = []
        for j in range(num_primes):
            if bits[j] == 1:
                value *= primes[j]
                equation_str.append(str(primes[j]))
        
        numbers_in_set.add(value)
        
        if not equation_str:
            print(f"For m={i:2d} (binary {i:04b}), value = 1")
        else:
            print(f"For m={i:2d} (binary {i:04b}), value = {' * '.join(equation_str)} = {value}")

    # The value of N is the maximum number in the set
    N = max(numbers_in_set)
    
    print("\nThe set of 16 numbers is:", sorted(list(numbers_in_set)))
    print("\nThe largest number in this set is N.")
    print(f"\nFinal Answer: N = {N}")

calculate_n()
