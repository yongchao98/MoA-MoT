import math

def solve_ratio():
    """
    Calculates the maximum achievable ratio of bidirectional conical power
    to the intensity along a line for the given physical system.
    """
    # The derived formula for the maximum ratio is:
    # Ratio = 4 * pi * (4/3 - 7 / (6 * sqrt(2)))

    # Define the constants from the formula
    c1 = 4.0
    pi = math.pi
    term1_num = 4.0
    term1_den = 3.0
    term2_num = 7.0
    term2_den_c1 = 6.0
    term2_den_c2 = 2.0
    
    # Calculate the components of the formula
    term1 = term1_num / term1_den
    sqrt_of_2 = math.sqrt(term2_den_c2)
    term2 = term2_num / (term2_den_c1 * sqrt_of_2)
    
    # Calculate the final ratio
    parenthesis_value = term1 - term2
    ratio = c1 * pi * parenthesis_value

    # Print the equation with its numerical components
    print("The final equation for the maximum ratio is:")
    print(f"{int(c1)} * pi * ({int(term1_num)}/{int(term1_den)} - {int(term2_num)} / ({int(term2_den_c1)} * sqrt({int(term2_den_c2)})))")
    print("\nSubstituting the values of the constants:")
    print(f"pi = {pi}")
    print(f"4 = {4}")
    print(f"3 = {3}")
    print(f"7 = {7}")
    print(f"6 = {6}")
    print(f"sqrt(2) = {sqrt_of_2}")
    
    # Print the final result
    print("\nThe calculated maximum ratio is:")
    print(ratio)

solve_ratio()