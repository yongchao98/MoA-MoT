import math

def calculate_average_constellation_size():
    """
    Calculates the theoretical average number of stars per constellation
    based on the nearest neighbor graph model.

    The formula is derived to be 8 / pi.
    """
    
    numerator = 8
    denominator = math.pi
    
    average_size = numerator / denominator
    
    print("The mathematical model of stars connected to their nearest neighbor leads to a theoretical average constellation size.")
    print("The derived formula for this average size is 8 / pi.")
    print("\nHere is the calculation:")
    print(f"{numerator} / {denominator} = {average_size}")

if __name__ == "__main__":
    calculate_average_constellation_size()