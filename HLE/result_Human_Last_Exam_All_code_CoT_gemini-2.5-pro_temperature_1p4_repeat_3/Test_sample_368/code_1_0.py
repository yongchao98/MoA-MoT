import math

def calculate_probability():
    """
    Calculates the probability that the given geometric series converges by iterating
    through all possible values of X, Y, and Z.
    """
    # The number of choices for each variable.
    # X can be any integer from -9 to 9, excluding 0.
    x_range = list(range(-9, 0)) + list(range(1, 10))
    # Y and Z can be any integer from 0 to 9.
    y_range = range(10)
    z_range = range(10)

    # Calculate the total number of possible combinations.
    total_count = len(x_range) * len(y_range) * len(z_range)

    # Initialize a counter for the combinations that lead to convergence.
    convergent_count = 0

    # The condition for convergence is |r| < 1, where r = 20k^2 + 24k.
    # This holds true for k in the union of two intervals.
    # We calculate the precise boundaries of these intervals using math.sqrt.
    lower_bound_1 = (-6 - math.sqrt(41)) / 10
    upper_bound_1 = (-6 - math.sqrt(31)) / 10
    lower_bound_2 = (-6 + math.sqrt(31)) / 10
    upper_bound_2 = (-6 + math.sqrt(41)) / 10

    # Iterate through every possible combination of X, Y, and Z.
    for x in x_range:
        for y in y_range:
            for z in z_range:
                # Calculate k for the current combination.
                # k = X + (10Y + 11Z) / 100
                k = x + (10 * y + 11 * z) / 100
                
                # Check if the calculated k falls within either of the two
                # intervals for convergence.
                is_in_interval1 = (lower_bound_1 < k < upper_bound_1)
                is_in_interval2 = (lower_bound_2 < k < upper_bound_2)

                if is_in_interval1 or is_in_interval2:
                    convergent_count += 1
    
    # The problem asks to output the final equation for the probability.
    print("The final equation for the probability is:")
    print(f"Probability = {convergent_count} / {len(x_range)} * {len(y_range)} * {len(z_range)}")
    print(f"Probability = {convergent_count} / {total_count}")

# Run the calculation and print the result.
calculate_probability()