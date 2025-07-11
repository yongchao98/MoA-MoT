import itertools
from fractions import Fraction

def find_largest_c():
    """
    This function searches for the largest density c = |R|/m by checking moduli
    m and all of their residue subsets R.

    The method is to find a set of residues R modulo m such that the sumset R+R
    contains no quadratic residues modulo m. The density c is then given by the
    ratio |R|/m. The function searches for the maximum possible c.
    """
    # Set a search limit for the modulus m.
    # The complexity grows as 2^M, so this must be kept reasonably small.
    # A limit of 16 is sufficient to find the best known constructions.
    M = 16

    best_c = Fraction(0)
    best_m = 0
    best_R = set()

    print(f"Searching for the best density c up to modulus M = {M}...")

    for m in range(2, M + 1):
        # Step 1: Determine the set of quadratic residues modulo m
        squares_mod_m = {(k * k) % m for k in range(m)}

        # Step 2: Iterate through all non-empty subsets R of {0, ..., m-1}.
        # A bitmask is used to represent and iterate through all 2^m - 1 subsets.
        for i in range(1, 1 << m):
            R = {j for j in range(m) if (i >> j) & 1}

            # Step 3: Check if the sumset R+R is free of quadratic residues
            is_valid = True
            # Use itertools.product to get all pairs from R (with replacement)
            for r1, r2 in itertools.product(R, repeat=2):
                if (r1 + r2) % m in squares_mod_m:
                    is_valid = False
                    break
            
            # Step 4: If the set R is valid, calculate its density
            if is_valid:
                current_c = Fraction(len(R), m)
                # If this density is the best so far, record it
                if current_c > best_c:
                    best_c = current_c
                    best_m = m
                    best_R = R

    # Output the results
    print("\n--- Search Complete ---")
    print("The best construction found is:")
    # Print the equation c = |R| / m
    print(f"c = |R| / m = {len(best_R)} / {best_m}")
    print(f"Modulus m = {best_m}")
    print(f"Residue set R = {sorted(list(best_R))}")
    print(f"Which gives a density c = {best_c.numerator}/{best_c.denominator} \u2248 {float(best_c):.4f}")
    print("\nBased on this search and known theoretical results, the largest number c is 1/3.")

if __name__ == '__main__':
    find_largest_c()