import math

def solve_it3_formulation():
    """
    Determines and prints the mathematical formulation for the vertical
    cross-section of a Gaussian-based Interval Type-3 Membership Function (IT3 MF).
    """

    # --- Define Parameters ---
    # These are example values to create a concrete formulation.

    # Parameters for the primary Gaussian membership function (upper bound of the main FOU)
    c = 5.0       # Center of the primary Gaussian
    sigma_x = 2.0 # Standard deviation of the primary Gaussian (along the x-axis)

    # Parameters for the vertical slice (secondary membership function)
    alpha = 0.9   # Scaling factor (height) for the vertical slice Gaussian
    sigma_u = 0.1 # Standard deviation for the vertical slice Gaussian (blurriness)

    # The fixed value of the primary input variable x where the vertical slice is taken
    x_prime = 3.0

    # --- Step 1: Calculate the Mean of the Vertical Slice Gaussian ---
    # The mean is the value of the primary upper membership function at x = x'.
    # The formula is: μ_upper(x') = exp( -0.5 * ( (x' - c) / σ_x )^2 )
    mu_upper_at_x_prime = math.exp(-0.5 * (((x_prime - c) / sigma_x) ** 2))

    # --- Step 2: Construct and Print the Final Formulation ---
    print("In the framework of an Interval Type-3 Fuzzy System, a vertical cross-section at a fixed input x=x' is a Type-2 Fuzzy Set.")
    print("Its membership bounds are often modeled using Gaussian functions.")
    print("\nThe formulation for the upper bound of a vertical cross-section, f(x', u), is as follows:")

    # Print the general symbolic formula for clarity
    print("\nGeneral Formula:")
    print("f(x', u) = α * exp( -0.5 * ( (u - μ_upper(x')) / σ_u )^2 )")
    print("  where μ_upper(x') = exp( -0.5 * ( (x' - c) / σ_x )^2 )")


    # Print the final, complete equation with the specific numbers calculated above.
    # This fulfills the requirement to output each number in the final equation.
    print("\n--- Final Equation with Example Values ---")
    print(f"For x' = {x_prime}, the calculated primary membership μ_upper({x_prime}) is {mu_upper_at_x_prime:.4f}.")
    print("\nThis results in the following specific formulation for the vertical slice:")
    final_equation = f"f({x_prime}, u) = {alpha} * exp( -0.5 * ( (u - {mu_upper_at_x_prime:.4f}) / {sigma_u} )^2 )"
    print(final_equation)

    print("\nWhere:")
    print("  - u: Represents the secondary input variable (the primary membership value).")
    print(f"  - α = {alpha}: The scaling factor determining the height of the vertical slice.")
    print(f"  - μ_upper({x_prime}) = {mu_upper_at_x_prime:.4f}: The mean of the vertical Gaussian, corresponding to the primary membership at x'.")
    print(f"  - σ_u = {sigma_u}: The standard deviation representing the uncertainty (blurriness) in the vertical dimension.")

# Execute the function to print the solution
solve_it3_formulation()
<<<f(3.0, u) = 0.9 * exp( -0.5 * ( (u - 0.6065) / 0.1 )^2 )>>>