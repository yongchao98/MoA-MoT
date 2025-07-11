import math

def solve_star_angle_problem():
    """
    Solves the relativistic star angle problem based on the provided information.

    The logic is derived from the properties of relativistic aberration. The target expression is:
    Ratio = (1 - cos(theta_14')) / (1 - cos(theta_34'))

    Through analysis of the Lorentz-invariant properties of light propagation, this ratio
    is shown to be equal to:
    Ratio = (1 - cos(theta_12')) / (1 - cos(theta_13'))
    where theta_12' and theta_13' are the angles given in the second frame of reference.
    """
    
    # Angles in the second frame of reference (in radians)
    theta_12_prime = math.pi / 2
    theta_13_prime = 3 * math.pi / 4

    # The numerator of the final expression is based on theta_12'
    cos_12_prime = math.cos(theta_12_prime)
    numerator_val = 1 - cos_12_prime

    # The denominator of the final expression is based on theta_13'
    cos_13_prime = math.cos(theta_13_prime)
    denominator_val = 1 - cos_13_prime
    
    # The final value is the ratio of these two quantities
    result = numerator_val / denominator_val
    
    print("The problem asks for the value of (1 - cos(theta_14')) / (1 - cos(theta_34')).")
    print("This can be shown to be equal to (1 - cos(theta_12')) / (1 - cos(theta_13')).\n")
    
    print("Let's calculate the numerator: (1 - cos(theta_12'))")
    print(f"The angle theta_12' between S1 and S2 is pi/2 radians.")
    print(f"cos(pi/2) = {cos_12_prime}")
    print(f"So, the numerator is 1 - {cos_12_prime:.1f} = {numerator_val:.1f}")
    print("-" * 30)

    print("Let's calculate the denominator: (1 - cos(theta_13'))")
    print(f"The angle theta_13' between S1 and S3 is 3*pi/4 radians.")
    print(f"cos(3*pi/4) = {cos_13_prime:.6f}")
    print(f"So, the denominator is 1 - ({cos_13_prime:.6f}) = {denominator_val:.6f}")
    print("-" * 30)
    
    print("The final equation is:")
    print(f"(1 - cos(theta_14')) / (1 - cos(theta_34')) = ({numerator_val:.1f}) / ({denominator_val:.6f})")
    print(f"Final Value = {result:.6f}")

solve_star_angle_problem()