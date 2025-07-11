def print_ac_loss_formula():
    """
    This function prints the normalized AC loss formula for an elliptical superconductor.

    The formula gives the loss per cycle per unit length (Q) as a function of
    the normalized current i = Im / Ic, where Im is the current amplitude and
    Ic is the critical current. The formula is valid for i < 1.

    The standard normalization for AC loss is 2*pi*Q / (mu_0 * Ic^2).
    """

    # The formula for the normalized loss, q(i)
    # This result is independent of the ellipse aspect ratio a/b for i < 1.
    formula = "2*pi*Q / (mu_0 * Ic^2) = (1-i)*ln(1-i) + (1+i)*ln(1+i) - i**2"

    print("The formula for the normalized AC loss per cycle for i < 1 is:")
    print(formula)

if __name__ == "__main__":
    print_ac_loss_formula()