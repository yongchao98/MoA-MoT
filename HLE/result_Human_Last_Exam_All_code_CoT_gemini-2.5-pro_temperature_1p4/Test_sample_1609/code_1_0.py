import math

def combine_terms(term1, term2):
    """
    Tries to combine two terms.
    A term is represented as a tuple (value, mask).
    e.g., minterm 5 (0101) is (5, 0).
    The term x-01 (represents xz'w) is (1, 8), i.e., (value=0001, mask=1000).
    Returns the new combined term, or None if they can't be combined.
    """
    # Masks must be the same to combine
    if term1[1] != term2[1]:
        return None
    
    # Values must differ by exactly one bit (a power of 2)
    diff = term1[0] ^ term2[0]
    if (diff & (diff - 1)) != 0:
        return None
        
    # Create the new term
    new_mask = term1[1] | diff
    new_value = term1[0] & term2[0]
    return (new_value, new_mask)

def count_prime_implicants(minterms, num_vars):
    """
    Calculates the number of prime implicants for a given set of minterms
    using the tabulation part of the Quine-McCluskey algorithm.
    """
    if not minterms:
        return 0
    # A function that is all 1s has one prime implicant: the constant "1"
    if len(minterms) == 2**num_vars:
        return 1
        
    # Start with the initial minterms as the first set of implicants
    current_terms = set((m, 0) for m in minterms)
    prime_implicants = set()
    
    while True:
        next_terms = set()
        merged_in_pass = set()
        
        # Iterate through all unique pairs of terms in the current set
        term_list = list(current_terms)
        for i in range(len(term_list)):
            for j in range(i + 1, len(term_list)):
                term1 = term_list[i]
                term2 = term_list[j]
                
                combined_term = combine_terms(term1, term2)
                if combined_term:
                    next_terms.add(combined_term)
                    merged_in_pass.add(term1)
                    merged_in_pass.add(term2)
        
        # Terms from the current set that were not merged are prime implicants
        unmerged = current_terms - merged_in_pass
        prime_implicants.update(unmerged)
        
        # If no new terms were created, the process is done
        if not next_terms:
            break
        
        # The new terms become the set for the next iteration
        current_terms = next_terms

    return len(prime_implicants)

def calculate_a4():
    """
    Finds a(4) by checking all self-dual functions of 4 variables.
    """
    n = 4
    num_total_minterms = 2**n
    num_pairs = num_total_minterms // 2
    
    # Create pairs of minterms for self-duality, e.g., (0, 15), (1, 14), ...
    pairs = [(i, num_total_minterms - 1 - i) for i in range(num_pairs)]
    
    max_pi_count = 0
    
    # A self-dual function is defined by choosing one minterm from each pair.
    # We can iterate through all 2^8 = 256 combinations.
    num_self_dual_functions = 2**num_pairs
    
    for i in range(num_self_dual_functions):
        minterms = set()
        
        # Build the set of minterms based on the i-th combination
        for k in range(num_pairs):
            if (i >> k) & 1:
                # Choose the second element of the pair
                minterms.add(pairs[k][1])
            else:
                # Choose the first element of the pair
                minterms.add(pairs[k][0])
                
        num_pis = count_prime_implicants(minterms, n)
        
        if num_pis > max_pi_count:
            max_pi_count = num_pis
            
    return max_pi_count

if __name__ == "__main__":
    a4 = calculate_a4()
    # The final value is the only "number in the final equation"
    print(a4)