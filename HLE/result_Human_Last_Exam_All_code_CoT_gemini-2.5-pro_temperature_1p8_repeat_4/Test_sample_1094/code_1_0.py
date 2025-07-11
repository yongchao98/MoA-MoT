# Plan:
# 1. State the physical problem: AC loss in an elliptic superconductor under a transport current.
# 2. Present the final formula for the normalized loss as a function of the normalized current `i`.
# 3. This formula, derived by W. T. Norris, is presented in its standard form.
# 4. The numerical constants '2' and '1' will be explicitly printed as part of the equation.
# 5. Explain all variables and constants in the formula.
# 6. Highlight the key insight that the result is independent of the ellipse's aspect ratio.

import math

# --- Introduction ---
print("This script provides the formula for the AC hysteresis loss in a superconducting bar of elliptic cross-section.")
print("The bar is in the critical state, characterized by a constant critical current density, Jc.")
print("The loss is calculated for an applied AC transport current with an amplitude Im less than the critical current Ic.")
print("-" * 80)

# --- Formula Presentation ---
# Define parts of the equation using Unicode for readability
normalized_loss_str = "2\u03C0Q/(\u03BC\u2080 I\u209C\u00B2)"
normalized_current_str = "i = I\u2098/I\u209C"

print(f"The normalized loss per cycle per unit length, {normalized_loss_str}, is given as a function of the normalized current amplitude, {normalized_current_str}, for i < 1.")
print("The formula is:")
print()

# The final equation. We print it with all its numerical components.
num_1 = 2
num_2 = 2
num_3 = 1
num_4 = 1

# Note: We are printing the string representation of the final formula, which inherently includes the numbers.
# 2*pi*Q/(mu_0*Ic^2) = (2 - i)*i + 2*(1 - i)*ln(1 - i)
final_equation = f"{normalized_loss_str} = ({num_1} - i)*i + {num_2}*({num_3} - i)*ln({num_4} - i)"
print(f"    {final_equation}")
print()
print("-" * 80)


# --- Variable Definitions ---
print("Where:")
print(f"  Q     : Hysteresis loss per cycle per unit length (in Joules/meter).")
print(f"  I\u2098    : Amplitude of the AC transport current.")
print(f"  I\u209C    : Critical current of the elliptic bar (I\u209C = Jc \u00B7 \u03C0 \u00B7 a \u00B7 b).")
print(f"  \u03BC\u2080   : Permeability of free space (4\u03C0 x 10\u207B\u2077 H/m).")
print(f"  ln()  : Natural logarithm.")
print(f"  a, b  : The semi-axes of the elliptical cross-section.")
print()

# --- Key Conclusion ---
print("Notably, this expression for the normalized AC loss is independent of the aspect ratio (a/b) of the elliptic bar.")
