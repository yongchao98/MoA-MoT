import math

def solve_earthquake_magnitude():
    """
    Calculates the expected maximum earthquake magnitude based on the derived formula.
    
    The formula is E = pi / (2 * log(2)) - 1, derived from the properties of
    Pareto and Logarithmic Series distributions.
    """
    
    # Define the constants from the final equation
    pi = math.pi
    two = 2.0
    natural_log_of_2 = math.log(2)
    one = 1.0

    # Print the explanation of the final formula and its components
    print("The problem reduces to calculating the value from the following equation:")
    print("E = pi / (2 * log(2)) - 1")
    print("\nWhere the numbers in the equation are:")
    print(f"pi = {pi}")
    print(f"2 = {two}")
    print(f"log(2) = {natural_log_of_2}")
    print(f"1 = {one}")

    # Calculate the result using the formula
    expected_max_magnitude = pi / (two * natural_log_of_2) - one

    # Print the final result
    print(f"\nThe expected maximum earthquake magnitude is: {expected_max_magnitude}")
    return expected_max_magnitude

if __name__ == "__main__":
    solve_earthquake_magnitude()