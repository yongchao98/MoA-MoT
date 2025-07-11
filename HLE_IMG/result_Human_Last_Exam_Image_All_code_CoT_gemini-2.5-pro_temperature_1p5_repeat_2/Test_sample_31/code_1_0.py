def calculate_rasmussen_invariant_7_4():
    """
    This script calculates the Rasmussen invariant for the knot shown in the image.
    The calculation is based on identifying the knot and applying a known theorem
    from knot theory.
    """

    # Step 1: Identify the knot
    knot_crossings = 7
    knot_index = 4
    knot_name = f"{knot_crossings}_{knot_index}"

    print("Step 1: Identify the knot from the diagram.")
    print(f"The knot in the image is the '{knot_name}' knot from the Alexander-Briggs notation.")
    print("-" * 30)

    # Step 2: State the relevant theorem
    print("Step 2: Use a theorem to simplify the calculation.")
    print(f"The {knot_name} knot is an 'alternating' knot.")
    print("A theorem by J. Rasmussen states that for any alternating knot K,")
    print("its Rasmussen s-invariant is equal to its signature, σ(K).")
    print("Formula: s(K) = σ(K)")
    print("-" * 30)

    # Step 3: Find the signature of the knot
    # This value is taken from standard knot theory tables.
    signature = -2
    print("Step 3: Find the signature of the knot.")
    print(f"From established knot theory tables, the signature of the {knot_name} knot is known.")
    print(f"σ({knot_name}) = {signature}")
    print("-" * 30)

    # Step 4: Calculate the Rasmussen invariant
    rasmussen_invariant = signature
    print("Step 4: Conclude the value of the Rasmussen invariant.")
    print(f"Applying the theorem from Step 2:")
    print(f"s({knot_name}) = σ({knot_name})")
    
    # Final equation with all numbers
    print("\nFinal Equation:")
    print(f"s({knot_crossings}_{knot_index}) = {rasmussen_invariant}")
    
    print("\nTherefore, the Rasmussen invariant of the knot is -2.")


if __name__ == "__main__":
    calculate_rasmussen_invariant_7_4()