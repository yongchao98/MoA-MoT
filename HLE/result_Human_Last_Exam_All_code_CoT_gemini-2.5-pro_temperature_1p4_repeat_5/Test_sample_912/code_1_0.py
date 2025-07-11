def solve_work_done():
    """
    This function provides a step-by-step derivation for the work done by the
    current source in the given magnetomechanical system over a cycle.
    """
    print("### Step-by-Step Derivation of Work Done by the Source ###")
    
    print("\n--- Step 1: Determine the Inductance L(x) ---")
    print("The magnetic circuit's reluctance changes with the block's position 'x'.")
    print("The circuit has two parallel magnetic paths within the main gap area:")
    print("  - Path 1 (through the block): Reluctance R_block = g / (μ * w * x)")
    print("  - Path 2 (through air): Reluctance R_air = g / (μ_0 * w * (D - x))")
    print("The total reluctance R_total(x) is the parallel combination:")
    print("  1 / R_total(x) = 1 / R_block + 1 / R_air = (μ*w*x)/g + (μ_0*w*(D-x))/g")
    print("  1 / R_total(x) = (w/g) * [ (μ - μ_0)x + μ_0*D ]")
    print("The inductance L(x) is given by L(x) = N^2 / R_total(x):")
    print("  L(x) = (N^2 * w / g) * [ (μ - μ_0)x + μ_0*D ]")

    print("\n--- Step 2: Calculate the Work for Each Path of the Cycle ---")
    print("The work done by the source is W = ∮ I dλ, where λ = L(x)I.")
    print("We evaluate the work for each of the four steps in the I-x plane:")
    print("  a) x from x1 to x2 (I = I1): W1 = ∫ I1 * d(L(x)I1) = I1^2 * (L(x2) - L(x1))")
    print("  b) I from I1 to I2 (x = x2): W2 = ∫ I*L(x2)*dI = (1/2) * L(x2) * (I2^2 - I1^2)")
    print("  c) x from x2 to x1 (I = I2): W3 = ∫ I2 * d(L(x)I2) = I2^2 * (L(x1) - L(x2))")
    print("  d) I from I2 to I1 (x = x1): W4 = ∫ I*L(x1)*dI = (1/2) * L(x1) * (I1^2 - I2^2)")

    print("\n--- Step 3: Sum the Work from All Steps ---")
    print("Total work W_cycle = W1 + W2 + W3 + W4.")
    print("Combining terms:")
    print("  W_cycle = [I1^2*(L(x2)-L(x1))] - [I2^2*(L(x2)-L(x1))] + [(1/2)*L(x2)*(I2^2-I1^2)] - [(1/2)*L(x1)*(I2^2-I1^2)]")
    print("  W_cycle = (L(x2)-L(x1)) * (I1^2-I2^2) + (1/2)*(L(x2)-L(x1)) * (I2^2-I1^2)")
    print("  W_cycle = -(L(x2)-L(x1)) * (I2^2-I1^2) + (1/2)*(L(x2)-L(x1)) * (I2^2-I1^2)")
    print("This simplifies to:")
    print("  W_cycle = -(1/2) * (L(x2) - L(x1)) * (I2^2 - I1^2)")

    print("\n--- Step 4: Substitute the Expression for L(x) ---")
    print("First, find the change in inductance, L(x2) - L(x1):")
    print("  L(x2) - L(x1) = (N^2*w/g) * (μ - μ_0) * (x2 - x1)")
    print("Now, substitute this into the equation for W_cycle:")
    print("  W_cycle = -(1/2) * [(N^2*w/g) * (μ - μ_0) * (x2 - x1)] * (I2^2 - I1^2)")

    print("\n### Final Equation ###")
    print("The final expression for the work done by the current source per cycle is:")
    
    # Define and print the components of the final formula
    numerator = "N^2 * w * (μ - μ_0) * (x2 - x1) * (I2^2 - I1^2)"
    denominator = "2g"
    print(f"W = - ( {numerator} ) / ( {denominator} )")

    # Re-arranging to match the format in the options
    print("\nOr, rearranged to match the options:")
    part1 = "(μ - μ_0)"
    part2 = "2g"
    part3 = "N^2 * w"
    part4 = "(x2 - x1)"
    part5 = "(I2^2 - I1^2)"
    print(f"W = - ( {part1} / {part2} ) * {part3} * {part4} * {part5}")


if __name__ == "__main__":
    solve_work_done()
