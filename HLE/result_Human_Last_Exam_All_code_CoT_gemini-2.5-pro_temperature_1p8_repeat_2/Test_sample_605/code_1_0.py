import sys

def calculate_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordstrom invariant for a Calabi-Yau Link.
    """
    # Weights for the variables (z1, z2, z3, z4, z5)
    weights = [22, 29, 49, 50, 75]
    n = len(weights)

    # For a Calabi-Yau manifold, the degree d is the sum of the weights.
    # This is a defining condition that resolves the apparent typo in the polynomial.
    d = sum(weights)

    print("To find the Crawley-Nordström invariant, we use the formula:")
    print("c = Σ (1 - 2*w_i / d) for i=1 to n\n")

    print(f"Given weights w = {tuple(weights)}")
    print(f"The number of variables is n = {n}")
    
    # Show the calculation of the degree d
    sum_of_weights_str = " + ".join(map(str, weights))
    print(f"For a Calabi-Yau link, the degree d equals the sum of the weights:")
    print(f"d = {sum_of_weights_str} = {d}\n")
    
    # --- Step 1: Substitute values into the formula ---
    step1_terms = [f"(1 - 2*{w}/{d})" for w in weights]
    print("Step 1: Substitute the values into the formula")
    print("c = " + " + ".join(step1_terms))
    print("")

    # --- Step 2: Simplify the 2*w_i part ---
    step2_numerators = [2 * w for w in weights]
    step2_terms = [f"(1 - {num}/{d})" for num in step2_numerators]
    print("Step 2: Perform the multiplication in the numerator")
    print("c = " + " + ".join(step2_terms))
    print("")

    # --- Step 3: Combine terms inside the parentheses ---
    step3_numerators = [d - num for num in step2_numerators]
    step3_terms = [f"({num}/{d})" for num in step3_numerators]
    print("Step 3: Perform the subtraction inside the parentheses")
    print("c = " + " + ".join(step3_terms))
    print("")

    # --- Step 4: Combine into a single fraction ---
    step4_numerator_str = " + ".join(map(str, step3_numerators))
    print("Step 4: Combine into a single fraction")
    print(f"c = ({step4_numerator_str}) / {d}")
    print("")

    # --- Step 5: Sum the numerators ---
    sum_of_numerators = sum(step3_numerators)
    print("Step 5: Sum the numbers in the numerator")
    print(f"c = {sum_of_numerators} / {d}")
    print("")
    
    # --- Step 6: Final Division ---
    final_result = sum_of_numerators / d
    print("Step 6: Perform the final division")
    # Use sys.stdout.write to avoid adding a newline, matching the requested format.
    sys.stdout.write(f"c = {int(final_result)}")


calculate_crawley_nordstrom_invariant()
# The final answer is enclosed below as requested.
# The code above calculates and prints the derivation.
print("\n") # Add a newline for separation before the final answer tag
print("<<<3>>>")