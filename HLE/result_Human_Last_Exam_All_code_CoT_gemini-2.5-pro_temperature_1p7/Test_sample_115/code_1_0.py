import math

def calculate_alpha():
    """
    Calculates the value of alpha, which is the golden ratio.
    """
    phi = (1 + math.sqrt(5)) / 2
    
    # The final equation is simply alpha = (1 + sqrt(5))/2
    print("The value of alpha is derived from the recurrence relation c_n \\approx \phi * c_{n-1}.")
    print("The base of the exponential growth alpha is the golden ratio phi.")
    print(f"alpha = (1 + sqrt(5)) / 2 = {phi}")

calculate_alpha()
