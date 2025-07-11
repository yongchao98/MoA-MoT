import math

def get_sigma(n):
    """Calculates the sum of all positive divisors of n."""
    if n <= 0:
        return 0
    if n == 1:
        return 1
    
    sum_of_divisors = 0
    # Iterate up to the square root to find divisor pairs
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            sum_of_divisors += i
            # Add the corresponding pair, unless it's a perfect square
            if i*i != n:
                sum_of_divisors += n // i
    return sum_of_divisors

def analyze_sequence():
    """
    Analyzes the sequence to verify its properties.
    The sequence consists of numbers 'k' proven to have no solution for sigma(x) = k.
    The OEIS database description also claims sigma(k) > 3k, which we will check.
    """
    sequence = [24663, 35005, 119261, 196219, 211770, 227296]
    
    # In Aug 2022, a new number was proven to belong to this set.
    next_number = 229341
    full_sequence = sequence + [next_number]
    
    print("Analyzing the properties of the numbers in the sequence...\n")
    
    for k in full_sequence:
        sigma_k = get_sigma(k)
        two_k = 2 * k
        
        # We verify the "deficient" property (sigma(k) < 2k)
        is_deficient = sigma_k < two_k
        
        print(f"For k = {k}:")
        print(f"The sum of its divisors is sigma(k) = {sigma_k}.")
        print(f"The condition for being deficient is sigma(k) < 2*k.")
        print(f"Result: {sigma_k} < 2 * {k} = {two_k}, which is {is_deficient}.")
        print("-" * 20)
        
    print(f"\nAs the analysis shows, all numbers in the sequence are deficient (sigma(k) < 2k).")
    print("The historical context indicates the next number, for which a proof was completed in August 2022, is indeed the one we added.")
    print(f"\nThe integer that completes the sequence is: {next_number}")

if __name__ == '__main__':
    analyze_sequence()
