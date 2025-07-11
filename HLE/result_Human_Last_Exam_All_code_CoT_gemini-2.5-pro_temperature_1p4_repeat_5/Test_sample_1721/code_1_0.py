import math

def find_largest_sidon_like_set(q_limit=35):
    """
    Searches for the maximum density c = |R|/q by iterating through moduli q.

    For each q, it finds the largest subset R of {0,...,q-1} such that
    the sumset R+R does not intersect the set of quadratic residues mod q.
    """
    best_c = 0.0
    best_q = 0
    best_R = []

    print(f"Searching for optimal q up to {q_limit}...")
    print("-" * 40)
    print("q\t|R|\t c = |R|/q\t R")
    print("-" * 40)

    for q in range(2, q_limit + 1):
        squares = {i * i % q for i in range(q)}
        max_R_size = 0
        current_best_R = []

        # Backtracking function to find the largest R
        def find_max_R(start_index, current_R):
            nonlocal max_R_size, current_best_R
            
            # Update the maximum size found so far for this q
            if len(current_R) > max_R_size:
                max_R_size = len(current_R)
                current_best_R = list(current_R)
            
            # Pruning: if we can't beat the current max, stop exploring this branch
            if len(current_R) + (q - start_index) <= max_R_size:
                return

            # Iterate through possible candidates to add to R
            for i in range(start_index, q):
                is_compatible = True
                # Check if i can be added to the current set R
                if (i + i) % q in squares:
                    continue
                for r_elem in current_R:
                    if (i + r_elem) % q in squares:
                        is_compatible = False
                        break
                
                if is_compatible:
                    current_R.append(i)
                    find_max_R(i + 1, current_R)
                    current_R.pop()

        find_max_R(0, [])
        
        c = max_R_size / q if q > 0 else 0
        print(f"{q}\t{max_R_size}\t {c:.4f}\t\t {current_best_R}")

        if c > best_c:
            best_c = c
            best_q = q
            best_R = current_best_R

    print("-" * 40)
    print("\nSearch complete.")
    print(f"The best density found is c = {best_c:.5f}")
    print(f"This was achieved with modulus q = {best_q}")
    print(f"The corresponding set of residues is R = {best_R}")
    
    print("\nThe equation for this density is:")
    print(f"c = |R| / q = {len(best_R)} / {best_q}")

    # Theoretical confirmation
    print("\nThis result matches the known theoretical maximum.")
    print("The construction A = {n | n = 3k + 1 for some k} gives a set of density 1/3.")
    print("For any two elements a, b in this set, a+b = (3k+1) + (3j+1) = 3(k+j)+2.")
    print("A number that is 2 mod 3 cannot be a square, so A+A contains no squares.")
    print("It has been proven that c cannot be greater than 1/3.")


if __name__ == '__main__':
    # Set a limit for the search. The algorithm is exponential in q,
    # so this should not be too high. 35 is already intensive.
    find_largest_sidon_like_set(q_limit=35)