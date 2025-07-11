import sympy

def calculate_alpha_exponent():
    """
    Calculates the specific heat scaling exponent alpha (α) using the first-order
    epsilon expansion in renormalization group theory.
    """
    # Define the parameters
    d = 3  # Spatial dimension
    d_c = 4  # Upper critical dimension
    n = 1  # Number of order parameter components (for the Ising model)

    # Calculate epsilon
    epsilon = d_c - d

    # Define symbols for a clean print output
    n_sym, eps_sym = sympy.symbols('n epsilon')
    alpha_formula_sym = (4 - n_sym) / (2 * (n_sym + 8)) * eps_sym

    # Perform the calculation using sympy for fractional results
    alpha_value = sympy.Rational((4 - n), (2 * (n + 8))) * epsilon

    # Print the explanation and results
    print("This script calculates the specific heat scaling exponent α using the epsilon expansion.")
    print("-" * 70)
    print("Step 1: Define the parameters.")
    print(f"  - Spatial dimension, d = {d}")
    print(f"  - Upper critical dimension, d_c = {d_c}")
    print(f"  - Order parameter components (Ising model), n = {n}\n")

    print("Step 2: Calculate the expansion parameter epsilon (ϵ).")
    print(f"  - Formula: ϵ = d_c - d")
    print(f"  - Calculation: ϵ = {d_c} - {d} = {epsilon}\n")

    print("Step 3: Calculate the scaling exponent alpha (α) to first order in ϵ.")
    print(f"  - General Formula: α = {alpha_formula_sym}")
    print(f"  - Plugging in the values (n={n}, ϵ={epsilon}):")

    # This loop constructs and prints the final equation with numbers
    numerator = 4 - n
    denominator_part1 = 2
    denominator_part2 = n + 8
    print(f"    α = ({4} - {n}) / ({denominator_part1} * ({n} + {8})) * {epsilon}")
    print(f"    α = {numerator} / ({denominator_part1} * {denominator_part2})")
    print(f"    α = {numerator} / {denominator_part1 * denominator_part2}")
    print(f"    α = {alpha_value}\n")

    print("Final Result:")
    print(f"The scaling exponent α for d={d} and n={n} is {alpha_value}.")
    print(f"As a decimal, this is approximately: {float(alpha_value):.4f}")

# Execute the function
calculate_alpha_exponent()