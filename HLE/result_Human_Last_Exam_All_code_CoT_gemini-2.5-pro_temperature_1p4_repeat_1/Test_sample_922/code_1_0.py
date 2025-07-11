import sys

def get_divisors_sum(n):
    """Calculates the sum of all positive divisors of n."""
    if n <= 0:
        return 0
    if n == 1:
        return 1
    
    # Start with 1 and n as divisors
    total = 1 + n
    # Iterate from 2 up to the square root of n
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            total += i
            # Add the corresponding divisor if it's not the square root
            if i * i != n:
                total += n // i
    return total

def find_next_term_in_sequence(start_n):
    """
    Finds the next integer n >= start_n such that 
    sigma(n) = sigma(n+1).
    """
    n = start_n
    print(f"Searching for the next term starting from n = {n}...")
    print("(This may take a few minutes)")
    
    while True:
        # Provide progress feedback to the user, printing to stderr
        if n % 5000 == 0:
            print(f"Still searching... currently at n = {n}", file=sys.stderr)

        # Calculate the sum of divisors for n and n+1
        sigma_n = get_divisors_sum(n)
        sigma_n_plus_1 = get_divisors_sum(n + 1)
        
        # Check if they are equal
        if sigma_n == sigma_n_plus_1:
            print(f"\nFound the next term: {n}")
            print(f"Verification of the condition sigma(n) = sigma(n+1):")
            # Outputting each number in the final equation as requested
            print(f"Sum of divisors of {n} is {sigma_n}")
            print(f"Sum of divisors of {n + 1} is {sigma_n_plus_1}")
            return n

        n += 1

if __name__ == "__main__":
    # The last number in the sequence provided by the user
    last_term_in_sequence = 227296
    
    # Start the search from the number immediately following the last known term
    next_term = find_next_term_in_sequence(last_term_in_sequence + 1)
