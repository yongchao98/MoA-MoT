import itertools

def solve_manifold_puzzle():
    """
    Solves the manifold puzzle by finding the lexicographically least tuple.
    The reasoning is embedded in the print statements.
    """

    print("--- Solving the Manifold Puzzle ---")

    print("\nStep 1: The 'Full' Condition")
    print("A compact, orientable n-manifold is 'full' if and only if its Euler characteristic is zero (χ(M) = 0).")
    print("This is because the existence of a metric of any odd index k (e.g., k=1) is obstructed by a non-zero Euler characteristic.")

    print("\nStep 2: Characterizing M(a,b)")
    print("The Euler characteristic of a surface of genus g is χ(M(g)) = 2 - 2g.")
    print("The Euler characteristic of a product of manifolds is the product of their characteristics.")
    print("Thus, for M(a,b) = M(a) x M(b), its characteristic is χ(M(a,b)) = (2 - 2a) * (2 - 2b).")
    print("For M(a,b) to be *not* full, we require χ(M(a,b)) ≠ 0, which means a ≠ 1 and b ≠ 1.")

    print("\nStep 3: The Connect-Sum Equation")
    print("Let X be the connect-sum of ℓ 4-manifolds M(aᵢ, bᵢ).")
    print("The Euler characteristic is χ(X) = (Σ χ(M(aᵢ, bᵢ))) - 2(ℓ - 1).")
    print("For X to be full, we need χ(X) = 0. This gives the key equation: Σ χᵢ = 2(ℓ - 1).")

    print("\nStep 4: Finding the Minimal ℓ")
    print("Each χᵢ = (2 - 2aᵢ)(2 - 2bᵢ) = 4(1 - aᵢ)(1 - bᵢ) is a multiple of 4.")
    print("Therefore, 2(ℓ - 1) must be a multiple of 4, which implies ℓ - 1 must be even.")
    print("This means ℓ must be an odd integer.")
    print("We are given that each M(aᵢ, bᵢ) is not full, so χᵢ ≠ 0. This rules out ℓ=1 (which would require χ₁=0).")
    print("The minimal possible value for ℓ is 3.")

    print("\nStep 5: The Equation for ℓ = 3")
    print("With ℓ = 3, our equation becomes: χ₁ + χ₂ + χ₃ = 2(3 - 1) = 4.")
    print("We must find pairs (a, b) with a,b ∈ {0, 2, 3,...} that solve this equation and form the lexicographically smallest tuple.")
    
    print("\nStep 6: Searching for the Solution")
    
    def chi_M(a, b):
        return (2 - 2 * a) * (2 - 2 * b)

    # Generate a sorted list of candidate pairs (a,b). A small search space is sufficient.
    max_genus_to_check = 5 
    candidate_pairs = []
    for a in range(max_genus_to_check):
        if a == 1: continue
        for b in range(max_genus_to_check):
            if b == 1: continue
            candidate_pairs.append((a, b))

    # Perform a systematic search for the three pairs.
    # By iterating through the sorted candidate_pairs and enforcing p1 <= p2 <= p3,
    # the first solution found is guaranteed to be the lexicographically minimal one.
    solution_found = False
    for p1 in candidate_pairs:
        for p2 in candidate_pairs:
            if p2 < p1: continue
            for p3 in candidate_pairs:
                if p3 < p2: continue
                
                chi1 = chi_M(p1[0], p1[1])
                chi2 = chi_M(p2[0], p2[1])
                chi3 = chi_M(p3[0], p3[1])

                if chi1 + chi2 + chi3 == 4:
                    print(f"Found the minimal set of pairs: {p1}, {p2}, {p3}")
                    
                    print("\nVerifying the final equation:")
                    print(f"χ{p1} + χ{p2} + χ{p3} = 4")
                    print(f"{chi1} + {chi2} + {chi3} = 4")
                    
                    final_tuple = p1 + p2 + p3
                    print("\nThe lexicographically least tuple is formed by concatenating these pairs:")
                    print(f"({','.join(map(str, final_tuple))})")
                    solution_found = True
                    break
            if solution_found:
                break
        if solution_found:
            break

if __name__ == '__main__':
    solve_manifold_puzzle()