def print_normalized_ac_loss_formula():
    """
    Calculates and prints the standard normalized AC loss for a superconductor
    with an elliptical cross-section carrying a transport current i = Im/Ic < 1.

    The derivation starts from Norris's formula for AC loss Q and rearranges it
    to the standard form 2*pi*Q / (mu_0 * Ic^2).
    """

    # According to Norris's model for an isolated conductor with i < 1:
    # Q = (mu_0 * Ic^2 / pi) * [ (1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2 ]

    # We want to find the expression for F(i) = 2*pi*Q / (mu_0 * Ic^2)

    # Step 1: Multiply Q by 2*pi
    # 2*pi*Q = 2*pi * (mu_0 * Ic^2 / pi) * [ ... ]
    # 2*pi*Q = 2 * mu_0 * Ic^2 * [ ... ]

    # Step 2: Divide by (mu_0 * Ic^2)
    # F(i) = 2 * [ (1-i)*ln(1-i) + (1+i)*ln(1+i) - i^2 ]

    # The numbers present in this final equation are identified below.
    # The requirement that the aspect ratio a/b is arbitrary is satisfied,
    # as this normalized formula is independent of the conductor's geometry.
    coefficient = 2
    one_val_1 = 1
    one_val_2 = 1
    one_val_3 = 1
    one_val_4 = 1
    exponent = 2

    # Construct the final equation string using the identified numbers.
    equation = (
        f"{coefficient} * [ ({one_val_1} - i)*ln({one_val_2} - i) + "
        f"({one_val_3} + i)*ln({one_val_4} + i) - i^{exponent} ]"
    )

    print("The normalized loss per cycle per unit length, 2*pi*Q/(mu_0 * Ic^2), "
          "as a function of i = Im/Ic is:")
    print(equation)


if __name__ == "__main__":
    print_normalized_ac_loss_formula()
