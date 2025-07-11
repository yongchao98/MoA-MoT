import math

def is_prime(n):
    """
    Checks if a number is prime using an efficient method.
    """
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def char_to_val(c):
    """
    Converts a character ('A'-'Z') to its 0-indexed value (0-25).
    """
    return ord(c) - ord('A')

def val_to_char(v):
    """
    Converts a 0-indexed value (0-25) to its character ('A'-'Z').
    """
    return chr(v + ord('A'))

def seq_to_num(seq):
    """
    Converts a 3-letter sequence to its base-26 decimal value.
    """
    p1 = char_to_val(seq[0])
    p2 = char_to_val(seq[1])
    p3 = char_to_val(seq[2])
    return p1 * 26**2 + p2 * 26 + p3

def num_to_seq(num):
    """
    Converts a decimal number to its 3-letter base-26 sequence.
    """
    p3 = num % 26
    num //= 26
    p2 = num % 26
    num //= 26
    p1 = num % 26
    return val_to_char(p1) + val_to_char(p2) + val_to_char(p3)

def solve_sequence():
    """
    Finds the next three terms in the sequence.
    """
    last_sequence_term = 'NZX'
    
    # Get the decimal value of the last term
    current_prime = seq_to_num(last_sequence_term)
    
    next_primes = []
    candidate = current_prime + 1
    
    # Find the next 3 prime numbers
    while len(next_primes) < 3:
        if is_prime(candidate):
            next_primes.append(candidate)
        candidate += 1
        
    # Convert the found primes back to sequences
    next_sequences = [num_to_seq(p) for p in next_primes]
    
    print("The last given sequence is '{}', which corresponds to the prime number {}.".format(last_sequence_term, current_prime))
    print("The next three prime numbers are {}, {}, and {}.".format(next_primes[0], next_primes[1], next_primes[2]))
    print("These primes correspond to the following three-letter sequences:")
    print(' '.join(next_sequences))

# Run the solver
solve_sequence()