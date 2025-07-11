import sys

# Redirect print to a string buffer to control the final output format if needed,
# but for this task, direct printing is fine.

def solve_homology_cobordism():
    """
    Calculates how many elements of the homology cobordism group can be represented
    by integral surgery on knots with at most four crossings.
    """
    print("This program determines how many elements of the homology cobordism group (Θ³_H)")
    print("can be represented by integral surgery on knots with at most four crossings.")
    print("\n---")

    print("Step 1: Define knots with at most 4 crossings and their Arf invariants.")
    # The Arf invariant of a knot K is 0 if its Alexander polynomial at -1, Δ_K(-1),
    # is congruent to ±1 mod 8, and 1 otherwise.
    # 0₁: Δ(t)=1, Δ(-1)=1. Arf=0.
    # 3₁: Δ(t)=t²-t+1, Δ(-1)=3. Arf=1.
    # 4₁: Δ(t)=-t+3-t⁻¹, Δ(-1)=5. Arf=1.
    knots = {
        "0₁ (Unknot)": 0,
        "3₁ (Trefoil)": 1,
        "4₁ (Figure-eight)": 1
    }
    print("The prime knots with at most 4 crossings and their Arf invariants are:")
    for name, arf in knots.items():
        print(f"  - {name}: Arf Invariant = {arf}")

    print("\n---")
    print("Step 2: Calculate achievable Rohlin invariants from integral surgeries.")
    
    # This set will store the unique Rohlin invariants we can generate.
    represented_invariants = set()

    # Process the unknot as a special case
    knot_name = "0₁ (Unknot)"
    print(f"\nProcessing knot: {knot_name}")
    print("  Any integral surgery (+1 or -1) on the unknot yields the 3-sphere S³.")
    # The Rohlin invariant of S³ is 0.
    rohlin_invariant_s3 = 0
    represented_invariants.add(rohlin_invariant_s3)
    print(f"  => Achieved Rohlin invariant: {rohlin_invariant_s3}. Current set of invariants: {represented_invariants}")

    # Process the other knots
    for name, arf in knots.items():
        if name == "0₁ (Unknot)":
            continue
        
        print(f"\nProcessing knot: {name} (Arf = {arf})")
        
        # For +1 surgery, μ = Arf(K)
        mu_plus_1 = arf % 2
        print(f"  - For +1 surgery: μ = Arf(K) mod 2 = {arf} mod 2 = {mu_plus_1}")
        represented_invariants.add(mu_plus_1)
        
        # For -1 surgery, μ = Arf(K) + 1
        mu_minus_1 = (arf + 1) % 2
        print(f"  - For -1 surgery: μ = (Arf(K) + 1) mod 2 = ({arf} + 1) mod 2 = {mu_minus_1}")
        represented_invariants.add(mu_minus_1)
        
        print(f"  => Current set of achievable invariants: {represented_invariants}")

    print("\n---")
    print("Step 3: Conclude and count the number of representable elements.")
    print("The homology cobordism group Θ³_H has 2 elements, distinguished by Rohlin invariants 0 and 1.")
    
    final_count = len(represented_invariants)
    final_set_str = str(sorted(list(represented_invariants))).replace('[', '{').replace(']', '}')

    print(f"The set of all Rohlin invariants that can be generated is {final_set_str}.")
    print("Since both possible invariants (0 and 1) can be generated, both elements of the group can be represented.")

    # The final "equation" showing the numbers involved in the count.
    print(f"\nFinal Equation: count(elements in {final_set_str}) = {final_count}")

solve_homology_cobordism()