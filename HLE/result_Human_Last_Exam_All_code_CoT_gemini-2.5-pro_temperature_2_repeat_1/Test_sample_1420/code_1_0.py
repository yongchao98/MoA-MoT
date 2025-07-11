import fractions

def calculate_alpha_exponent():
    """
    Calculates the specific heat exponent alpha (α) for a d=3 system
    using a first-order epsilon expansion.
    """
    # Parameters for the calculation
    d = 3      # Spatial dimension
    d_c = 4    # Upper critical dimension
    n = 1      # Number of order parameter components (assuming the Ising model)

    # Explanation of the method
    print("Calculating the specific heat exponent alpha (α) for a d=3 system.")
    print("The calculation uses a first-order epsilon expansion near the upper critical dimension.")
    print(f"We assume the system is in the Ising universality class (n={n} component order parameter).")
    print("-" * 50)

    # Step 1: Calculate epsilon
    epsilon = d_c - d
    print("Step 1: Calculate the expansion parameter epsilon (ϵ)")
    print(f"ϵ = d_c - d = {d_c} - {d} = {epsilon}")
    print("")

    # Step 2: Calculate alpha using the first-order formula
    print("Step 2: Use the formula for α to first order in ϵ")
    print("α = (4 - n) / (2 * (n + 8)) * ϵ")
    print(f"Plugging in n={n} and ϵ={epsilon}:")
    
    # Show the calculation step-by-step
    numerator = 4 - n
    denominator_part1 = 2
    denominator_part2 = n + 8
    denominator = denominator_part1 * denominator_part2
    
    print(f"α = ({4} - {n}) / ({denominator_part1} * ({n} + {8})) * {epsilon}")
    print(f"α = {numerator} / ({denominator_part1} * {denominator_part2}) * {epsilon}")
    print(f"α = {numerator} / {denominator} * {epsilon}")

    # Calculate final result as a fraction and a float
    alpha_fraction = fractions.Fraction(numerator, denominator) * epsilon
    alpha_float = float(alpha_fraction)
    
    print(f"α = {alpha_fraction.numerator}/{alpha_fraction.denominator}")
    print(f"As a decimal, α ≈ {alpha_float}")
    print("-" * 50)

# Run the calculation
calculate_alpha_exponent()