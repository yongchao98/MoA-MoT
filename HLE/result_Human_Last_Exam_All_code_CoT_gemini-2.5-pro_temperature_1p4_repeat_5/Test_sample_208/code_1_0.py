import math

def calculate_packing_radius():
    """
    Calculates the radius of a circle tightly packing 14 smaller circles of radius 1.
    """
    # Define the parameters for the problem
    # r is the radius of the small circles
    r = 1.0
    # k is the number of circles in the outer layer for the optimal packing of 14 circles
    k = 8

    # Explain the context and the formula being used
    print("This problem is a case of 'circle packing in a circle'.")
    print("For 14 circles of radius 1, the optimal arrangement has an outer layer of 8 circles.")
    print("The radius 'R' of the containing circle is calculated with the formula:")
    print("R = r * (1 + 1 / sin(pi / k))\n")

    # Perform the calculation step-by-step and print the equation with numbers
    print("Plugging in the values r=1 and k=8:")
    
    # Step 1: Show the initial formula with numbers
    print(f"R = {r} * (1 + 1 / sin(pi / {k}))")

    # Step 2: Calculate intermediate values to show the process
    sin_val = math.sin(math.pi / k)
    print(f"R = {r} * (1 + 1 / {sin_val:.7f})")
    
    inv_sin_val = 1 / sin_val
    print(f"R = {r} * (1 + {inv_sin_val:.7f})")
    
    R = r * (1 + inv_sin_val)
    print(f"R = {R:.7f}\n")

    # Format the final result to 4 significant digits
    # The '.4g' format specifier is used for rounding to a certain number of significant figures.
    formatted_R = f"{R:.4g}"

    print(f"The radius of the large circle, rounded to 4 significant digits, is: {formatted_R}")

if __name__ == "__main__":
    calculate_packing_radius()