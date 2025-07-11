import sys

# Define a function to handle the calculation and printing
def calculate_specific_heat_exponent():
    """
    Calculates the specific heat scaling exponent alpha using the first-order
    epsilon expansion for the 3D Ising model.
    """
    # Step 1: Define parameters
    d = 3      # Spatial dimension
    d_c = 4    # Upper critical dimension for the O(n) model
    n = 1      # Number of spin components for the Ising model

    print(f"Calculating the specific heat exponent α for a system with d = {d}.")
    print(f"We use the epsilon expansion for the O(n) model, assuming the Ising universality class (n = {n}).")
    print(f"The upper critical dimension for this model is d_c = {d_c}.")
    print("-" * 30)

    # Step 2: Calculate epsilon
    epsilon = d_c - d
    print("Step 1: Calculate the expansion parameter ϵ = d_c - d.")
    print(f"ϵ = {d_c} - {d} = {epsilon}")
    print("-" * 30)

    # Step 3: Use the first-order formula for alpha
    print("Step 2: Apply the first-order expansion formula for α:")
    print("α = (4 - n) / (2 * (n + 8)) * ϵ")
    print("-" * 30)

    # Step 4: Substitute values and calculate
    print("Step 3: Substitute the values and compute α.")
    
    # Breaking down the formula for clear output
    numerator = 4 - n
    denominator_part1 = 2
    denominator_part2 = n + 8
    denominator = denominator_part1 * denominator_part2
    
    # The final calculation
    alpha = (numerator / denominator) * epsilon

    print(f"α = ({4} - {n}) / ({denominator_part1} * ({n} + {8})) * {epsilon}")
    print(f"α = {numerator} / ({denominator_part1} * {denominator_part2}) * {epsilon}")
    print(f"α = {numerator} / {denominator} * {epsilon}")
    
    # To display the result as a fraction
    # In Python 3, a simple way is using the fractions module
    try:
        from fractions import Fraction
        alpha_fraction = Fraction(numerator, denominator) * epsilon
        print(f"α = {alpha_fraction.numerator}/{alpha_fraction.denominator}")
    except ImportError:
        # Fallback for older python or minimal environments
        pass
        
    print(f"α ≈ {alpha}")
    
    # Store the final numerical answer for the "answer" block
    # This part will not be printed to the user console
    sys.stdout = open(os.devnull, 'w')
    print(f"<<<{alpha}>>>")
    sys.stdout = sys.__stdout__


# Execute the function
calculate_specific_heat_exponent()
