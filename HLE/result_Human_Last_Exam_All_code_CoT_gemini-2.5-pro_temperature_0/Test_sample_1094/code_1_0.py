def solve_ac_loss_formula():
    """
    This function provides the analytical formula for the normalized AC loss
    in an elliptical superconductor for a transport current i = Im/Ic < 1.
    """

    # The problem asks for the loss per cycle Q in the normalized form:
    # 2 * pi * Q / (mu_0 * Ic^2)
    #
    # The Norris formula for loss Q is:
    # Q = (mu_0 * Ic^2 / pi) * [(1-i)*ln(1-i) + (1+i)*ln(1+i) - i**2]
    #
    # To get the normalized form, we rearrange the formula:
    # 2 * pi * Q = 2 * mu_0 * Ic^2 * [(1-i)*ln(1-i) + (1+i)*ln(1+i) - i**2]
    # (2 * pi * Q) / (mu_0 * Ic^2) = 2 * [(1-i)*ln(1-i) + (1+i)*ln(1+i) - i**2]

    # We will represent the variable i = Im/Ic as a string 'i' for the output.
    i = "i"

    # Construct the final equation string, showing each number and operation.
    final_equation = f"2 * ((1 - {i})*ln(1 - {i}) + (1 + {i})*ln(1 + {i}) - {i}**2)"

    print("The normalized AC loss per cycle per unit length, 2*pi*Q/(mu_0*Ic^2), as a function of i = Im/Ic is:")
    print(final_equation)

# Execute the function to print the result
solve_ac_loss_formula()