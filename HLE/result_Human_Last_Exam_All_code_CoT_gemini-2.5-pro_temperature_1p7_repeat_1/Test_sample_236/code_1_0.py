import sys

def solve():
    """
    Calculates the number of homology cobordism group elements representable
    by integral surgery on knots with at most four crossings.
    """
    
    print("Goal: Find how many elements of the homology cobordism group Θ^3_H = Z/2Z can be represented.")
    print("Method: Check representatives from integral surgeries on knots with <= 4 crossings.")
    print("An element is identified by its Rohlin invariant (0 for trivial, 1 for non-trivial).\n")
    print("Key Formulas:")
    print("1. Rohlin invariant of +1 surgery: μ(M(K, +1)) = 0 (always represents the trivial element).")
    print("2. Rohlin invariant of -1 surgery: μ(M(K, -1)) = Arf(K).")
    print("3. Arf(K) is 0 if Δ_K(-1) ≡ ±1 (mod 8), and 1 if Δ_K(-1) ≡ ±3 (mod 8).\n")
    
    # Knots with at most 4 crossings and their Alexander polynomials evaluated at t=-1.
    # We use Δ(-1) to calculate the Arf invariant.
    knots_info = {
        '0_1': {'name': 'Unknot', 'delta_neg_one': 1},
        '3_1': {'name': 'Trefoil Knot', 'delta_neg_one': -3},
        '4_1': {'name': 'Figure-Eight Knot', 'delta_neg_one': 5}
    }
    
    # This set will store the unique elements (0 or 1) that we can represent.
    represented_elements = set()
    
    # Element 0 (trivial) is always representable by +1 surgery on any knot.
    # We can use the unknot as a concrete example.
    represented_elements.add(0)
    print("--- Analyzing +1 Surgeries ---")
    print("Knot: 0_1 (Unknot)")
    print("Surgery: M(0_1, +1) results in a manifold with Rohlin invariant 0.")
    print("Equation: M(0_1, +1) => 0")
    print("This represents the trivial element '0'.\n")

    print("--- Analyzing -1 Surgeries ---")
    # Check -1 surgeries to see which elements they represent.
    for knot_id, info in knots_info.items():
        delta_val = info['delta_neg_one']
        
        # Calculate Arf invariant
        val_mod_8 = delta_val % 8
        if val_mod_8 in (1, 7): # ±1 mod 8
            arf_invariant = 0
        elif val_mod_8 in (3, 5): # ±3 mod 8
            arf_invariant = 1
        else:
            # This case should not be reached for knots.
            arf_invariant = -1 

        represented_elements.add(arf_invariant)
        
        print(f"Knot: {knot_id} ({info['name']})")
        print(f"Δ(-1) = {delta_val}, which is {val_mod_8} (mod 8).")
        print(f"The Arf invariant is {arf_invariant}.")
        print(f"Surgery: M({knot_id}, -1) represents the element '{arf_invariant}'.")
        print(f"Equation: M({knot_id}, -1) => {arf_invariant}")
        print("-" * 20)

    # Final summary
    print("\n--- Conclusion ---")
    print(f"The set of all representable elements is: {represented_elements}.")
    
    num_elements = len(represented_elements)
    
    print(f"Both the trivial element '0' and the non-trivial element '1' of the homology cobordism group can be represented.")
    print(f"Therefore, the total number of representable elements is {num_elements}.")

solve()
<<<2>>>