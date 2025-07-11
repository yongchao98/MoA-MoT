import math

def calculate_hitting_probability():
    """
    Calculates the probability that a 2D random walk visits a target set
    before leaving a large disk, using a continuous potential theory approximation.
    """
    # Parameters of the problem
    z0 = (0, 300)
    z1 = (0, 0)
    z2 = (2, 0)
    R = 1000.0
    r0 = 1.0  # Effective radius of a lattice point (cutoff)

    # Calculate the required distances
    d01 = math.sqrt((z0[0] - z1[0])**2 + (z0[1] - z1[1])**2)
    d02 = math.sqrt((z0[0] - z2[0])**2 + (z0[1] - z2[1])**2)
    d12 = math.sqrt((z1[0] - z2[0])**2 + (z1[1] - z2[1])**2)

    # Numerator of the probability formula
    # log( (|z0-z1|*|z0-z2|) / R^2 )
    numerator_arg = (d01 * d02) / (R**2)
    numerator = math.log(numerator_arg)

    # Denominator of the probability formula
    # log( (r0*|z1-z2|) / R^2 )
    denominator_arg = (r0 * d12) / (R**2)
    denominator = math.log(denominator_arg)
    
    # Calculate the probability
    probability = numerator / denominator

    # Output the equation with the calculated values
    print("The hitting probability P is calculated using the formula:")
    print("P = log( (d01 * d02) / R^2 ) / log( (r0 * d12) / R^2 )")
    print("\nWith the following values:")
    print(f"d01 (distance from start to target 1): {d01:.3f}")
    print(f"d02 (distance from start to target 2): {d02:.3f}")
    print(f"d12 (distance between targets): {d12:.3f}")
    print(f"R (disk radius): {R:.1f}")
    print(f"r0 (cutoff radius): {r0:.1f}")
    
    print("\nPlugging in the values:")
    print(f"P = log( ({d01:.3f} * {d02:.3f}) / {R**2:.1f} ) / log( ({r0:.1f} * {d12:.3f}) / {R**2:.1f} )")
    print(f"P = log({numerator_arg:.6f}) / log({denominator_arg:.8f})")
    print(f"P = {numerator:.4f} / {denominator:.4f}")

    print(f"\nThe final probability is: {probability:.3f}")

calculate_hitting_probability()