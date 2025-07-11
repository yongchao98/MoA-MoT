import sys

def find_probability_pm(m):
    """
    Calculates the probability P_m for a given positive integer m.

    Args:
      m: A positive integer.
    """
    if not isinstance(m, int) or m <= 0:
        print("Error: m must be a positive integer.", file=sys.stderr)
        return

    # Calculate the number of favorable pairs using the formula m(m+2)
    numerator_val = m * (m + 2)

    # Calculate the total number of pairs using the formula (2m+1)(4m+1)
    denominator_val = (2 * m + 1) * (4 * m + 1)

    # Calculate the final probability as a float
    probability_val = numerator_val / denominator_val

    # Print the results in a step-by-step equation format, showing each number involved.
    # This fulfills the request to "output each number in the final equation".
    print(f"The general formula for P_m is: m(m+2) / ((2m+1)(4m+1))")
    print(f"For m = {m}:")
    print(f"P_{m} = ({m} * ({m} + 2)) / ((2 * {m} + 1) * (4 * {m} + 1))")
    print(f"P_{m} = {numerator_val} / {denominator_val}")
    print(f"P_{m} = {probability_val}")


# --- Example Usage ---
# You can change the value of m here to calculate P_m for different values.
m_value = 1
find_probability_pm(m_value)

print("\n--- Another Example ---")
m_value = 2
find_probability_pm(m_value)