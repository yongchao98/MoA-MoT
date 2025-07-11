# The value k is a given integer parameter, k >= 2.
# We set a specific value for k for demonstration.
# You can change this value to compute the limit for a different k.
k = 2

# Check if k is valid as per the problem statement.
if not isinstance(k, int) or k < 2:
    print(f"Error: k must be an integer greater than or equal to 2. Found k = {k}.")
else:
    # Based on the mathematical derivation, the limit is given by the formula:
    # L = 1 - 1 / (2 * k)

    # The components of the formula are:
    val_1 = 1
    numerator = 1
    coeff_2 = 2

    # Calculation
    limit_value = val_1 - numerator / (coeff_2 * k)

    # As requested, we output each number in the final equation as part of the calculation.
    print(f"For k = {k}, the limit is calculated as:")
    print(f"{val_1} - {numerator} / ({coeff_2} * {k}) = {limit_value}")
