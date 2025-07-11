import math

def solve_integral():
    """
    This function calculates the integral of f(x,y,z) = z^2*(x^2+y^2) over a cone
    and prints the step-by-step derivation and final answer.
    """
    # --- Parameters ---
    H = 2  # Height of the cone
    R = 3  # Radius of the cone's base

    # --- Step-by-step Explanation ---
    print("### Solving the Integral of f(x,y,z) = z^2*(x^2+y^2) over a Cone ###\n")
    print("The integral is best solved in cylindrical coordinates (r, theta, z).")
    print(f"The cone has height H={H} and radius R={R}.")
    print("The integral is set up as: I = Integral(z=0 to 2) [ Integral(theta=0 to 2*pi) [ Integral(r=0 to 3/2*(2-z)) [ z^2 * r^3 dr ] d(theta) ] dz ]\n")

    print("Step 1: Integrate with respect to 'r'.")
    print("  Integral(z^2 * r^3)dr = z^2 * r^4 / 4")
    print("  Evaluating from r=0 to r=(3/2)*(2-z) gives: (81/64) * z^2 * (2-z)^4\n")

    print("Step 2: Integrate the result with respect to 'theta'.")
    print("  The expression has no 'theta', so we multiply by 2*pi.")
    print("  Result: 2*pi * (81/64) * z^2 * (2-z)^4 = (81*pi/32) * z^2 * (2-z)^4\n")

    print("Step 3: Integrate the result with respect to 'z'.")
    print("  We need to compute Integral from z=0 to 2 of [(81*pi/32) * z^2 * (2-z)^4] dz.")
    print("  The constant part is (81*pi/32).")
    # The value of the definite integral of z^2 * (2-z)^4 from 0 to 2 is 128/105.
    integral_z_part_num = 128
    integral_z_part_den = 105
    print(f"  The definite integral part evaluates to: {integral_z_part_num}/{integral_z_part_den}\n")

    print("Step 4: Calculate the final result.")
    # Components from the derivation
    coeff_num = 81
    coeff_den = 32

    # Final calculation
    final_num_unsimplified = coeff_num * integral_z_part_num
    final_den_unsimplified = coeff_den * integral_z_part_den
    
    # Simplify the fraction
    common_divisor = math.gcd(final_num_unsimplified, final_den_unsimplified)
    final_num_simplified = final_num_unsimplified // common_divisor
    final_den_simplified = final_den_unsimplified // common_divisor

    print("The final equation is the product of the constant part and the z-integral result:")
    print(f"  I = ( {coeff_num} * pi / {coeff_den} ) * ( {integral_z_part_num} / {integral_z_part_den} )")
    print(f"  I = ( {final_num_unsimplified} * pi ) / {final_den_unsimplified}")
    print(f"After simplifying the fraction, the final exact answer is:")
    print(f"  I = ({final_num_simplified} * pi) / {final_den_simplified}")
    
    # Calculate numerical value
    final_value = (final_num_simplified * math.pi) / final_den_simplified
    print(f"\nThe approximate numerical value of the integral is: {final_value}")

solve_integral()

# The final result in the requested format
final_numerator = 108
final_denominator = 35
final_answer_value = (final_numerator * math.pi) / final_denominator
# print(f'<<<{final_answer_value}>>>') # for internal check
# Let's print the fraction itself as the content as it's more precise
print(f'<<<(108*pi)/35>>>')