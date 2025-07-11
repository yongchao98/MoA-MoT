import math

def calculate_series_sum():
    """
    Calculates the value of the series sum based on the derived formula.
    The sum is S = 1 - (pi^2 / 12) + (ln(2)^2 / 2).
    """
    # Calculate each term of the formula
    pi_squared_over_12 = math.pi**2 / 12
    ln2_squared_over_2 = math.log(2)**2 / 2
    
    # Calculate the final sum
    result = 1 - pi_squared_over_12 + ln2_squared_over_2
    
    # Print the equation with the calculated values of each term
    print("The final equation is S = 1 - (pi^2 / 12) + (ln(2)^2 / 2)")
    print(f"Term 1: 1")
    print(f"Term 2 (pi^2 / 12): {pi_squared_over_12}")
    print(f"Term 3 ((ln(2))^2 / 2): {ln2_squared_over_2}")
    print(f"Final Value S = 1 - {pi_squared_over_12} + {ln2_squared_over_2} = {result}")

if __name__ == "__main__":
    calculate_series_sum()