import math

def solve_sum():
    """
    Calculates the value of the sum based on the derived formula.
    The sum is S = 1 - (pi^2 / 12) + (ln(2)^2 / 2).
    """
    
    # Calculate the components of the formula
    pi_squared_over_12 = math.pi**2 / 12
    ln2_squared_over_2 = math.log(2)**2 / 2
    
    # Calculate the final value
    final_value = 1 - pi_squared_over_12 + ln2_squared_over_2
    
    # Print the equation with the calculated numbers
    print("The symbolic form of the sum is: 1 - (pi^2 / 12) + (ln(2)^2 / 2)")
    print("\nCalculating the numerical value:")
    print(f"1 - {pi_squared_over_12} + {ln2_squared_over_2}")
    print(f"= {final_value}")

if __name__ == "__main__":
    solve_sum()