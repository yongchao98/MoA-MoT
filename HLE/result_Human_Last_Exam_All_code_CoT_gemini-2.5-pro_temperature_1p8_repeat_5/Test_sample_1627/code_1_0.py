import sys

def calculate_stress_ratio():
    """
    Calculates the ratio of maximum tangential stress to internal pressure
    for a thick-walled cylinder where the outer radius is twice the inner radius.
    """
    # The problem asks for the ratio of maximum tangential stress (sigma_t_max)
    # to internal pressure (P_i) for a thick-walled cylinder.
    #
    # The governing formula for this ratio is:
    # Ratio = (b^2 + a^2) / (b^2 - a^2)
    # where 'a' is the inner radius and 'b' is the outer radius.
    #
    # We are given that the outer radius is twice the inner radius, so:
    # b = 2 * a
    #
    # We can define a radius ratio k = b / a = 2.
    k = 2.0

    print("Step 1: Define the ratio of outer radius 'b' to inner radius 'a'.")
    print(f"k = b / a = {k}\n")

    print("Step 2: The ratio of maximum tangential stress to internal pressure is given by (k^2 + 1) / (k^2 - 1).")
    
    # Calculate the numerator and denominator of the final equation
    numerator = k**2 + 1
    denominator = k**2 - 1

    print("Step 3: Substitute k = 2 into the equation.")
    print(f"Numerator = k^2 + 1 = {k}^2 + 1 = {k**2} + 1 = {numerator}")
    print(f"Denominator = k^2 - 1 = {k}^2 - 1 = {k**2} - 1 = {denominator}\n")

    # The final ratio
    ratio_value = numerator / denominator

    print(f"The final equation for the ratio is: {int(numerator)} / {int(denominator)}")
    print(f"The numerical value of the ratio is: {ratio_value}")

# Execute the calculation
calculate_stress_ratio()

# Based on the calculation, the ratio is 5/3.
# Comparing this with the answer choices:
# A. 1.5:1  (1.5)
# B. 4:1    (4.0)
# C. 4.5:1  (4.5)
# D. 2:1    (2.0)
# E. 5/3    (1.666...)
# F. 5/9    (0.555...)
# G. 5:1    (5.0)
# H. 1.667  (Approximation of 5/3)
# I. 6:1    (6.0)
# J. 3:1    (3.0)
#
# The exact answer is 5/3, which corresponds to choice E.
# We will use sys.stdout.write to output the final choice without extra text from the print function.
# This part is just for providing the final answer in the specified format for the platform.
# The user will only see the output from the function above.
final_answer = 'E'
sys.stdout.write(f"\n<<<{final_answer}>>>")