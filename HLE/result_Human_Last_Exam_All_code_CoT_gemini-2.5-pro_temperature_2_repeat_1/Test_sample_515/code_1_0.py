import math

def potential_kernel(x, y):
    """
    Calculates the potential kernel a(x,y) using the asymptotic formula.
    a(v) approx 2/pi * log(||v||) + K
    K = (2*gamma + log(8))/pi
    We only need differences, or ratios of its value for large vectors,
    which makes the constant K less important for some calculations, but we include it for accuracy.
    """
    R = math.sqrt(x**2 + y**2)
    if R == 0:
        return 0.0
    
    # Constants
    gamma = 0.5772156649
    # K = (2 * gamma + math.log(8)) / math.pi
    
    # Let's compute a(v) = (2/pi)*log(||v||)+K and also a specific value a(1,0)=4/pi
    # which is used as a reference. Let's define a normalized potential where a(1)=1.
    # a_norm(x) = a(x)/a(1).
    # Then P = (1 + a_norm(z) - a_norm(x0-z))/2.
    # a(x)/a(1) = (2/pi*log(x)+K)/(4/pi) = 1/2*log(x)+pi*K/4
    
    # Let's use the formulas as derived.
    return (2 / math.pi) * math.log(R) + (2 * gamma + math.log(8)) / math.pi

def main():
    """
    Main function to calculate the required probability.
    """
    x0_x, x0_y = 0, 1
    z_x, z_y = 3600, 0

    # The formula for hitting probability of point z from x0 is:
    # P = (a(x0) + a(z) - a(x0 - z)) / (2 * a(x0))
    # where a(v) is the potential kernel.
    
    # For x0=(0,1), a(x0) = a(1) = 4/pi
    a_x0 = 4 / math.pi

    # For z=(3600,0)
    a_z = potential_kernel(z_x, z_y)
    
    # For x0-z=(-3600,1)
    a_x0_minus_z = potential_kernel(x0_x - z_x, x0_y - z_y)

    # Calculate the probability
    prob = (a_x0 + a_z - a_x0_minus_z) / (2 * a_x0)
    
    # Output the result formatted to two significant digits
    print(f"The starting point is x0 = ({x0_x}, {x0_y})")
    print(f"The target point is z = ({z_x}, {z_y})")
    print(f"The value of the potential kernel at x0=(0,1) is a(x0) = 4/\u03C0 \u2248 {a_x0:.4f}")
    print(f"The value of the potential kernel at z=(3600,0) is a(z) \u2248 {a_z:.4f}")
    print(f"The value of the potential kernel at x0-z=(-3600,1) is a(x0-z) \u2248 {a_x0_minus_z:.4f}")
    print(f"The probability is given by the formula P = (a(x0) + a(z) - a(x0-z))/(2*a(x0))")
    print(f"P \u2248 ({a_x0:.4f} + {a_z:.4f} - {a_x0_minus_z:.4f}) / (2 * {a_x0:.4f})")
    print(f"P \u2248 {prob:.4f}")
    
    # Final answer rounded to two significant digits.
    print(f"\nThe approximate probability is {prob:.2g}")


if __name__ == "__main__":
    main()
