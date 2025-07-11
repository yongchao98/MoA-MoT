import numpy as np

# In toroidal plasma physics, physical quantities that vary with the poloidal angle
# 'theta' are typically expanded in a Fourier series due to the periodic nature
# of this coordinate (0 to 2*pi).

# A function F(theta) can be represented as:
# F(theta) = A0 + sum_{m=1 to N} [ Am * cos(m*theta) + Bm * sin(m*theta) ]
# where 'm' is the poloidal mode number.

# Let's define some example coefficients for a hypothetical physical quantity.
# A0 represents the average value (the m=0 component).
# A1, B1 represent the m=1 "dipole" component.
# A2, B2 represent the m=2 "quadrupole" component.

coeff_A0 = 10.0
coeff_A1 = 2.5
coeff_B1 = 0.0
coeff_A2 = 0.0
coeff_B2 = -1.2

# The task is to show the final equation with each number printed out.
# This code will construct and print the mathematical expression for a
# function composed of m=0, m=1, and m=2 poloidal modes.

print("A physical quantity F(theta) in a toroidal system can be represented by a Fourier series.")
print("Here is an example equation for a function built from m=0, 1, and 2 modes:")
print("\nGeneral form: F(theta) = A0 + A1*cos(theta) + B1*sin(theta) + A2*cos(2*theta) + B2*sin(2*theta) + ...")

# Print the specific equation using our defined coefficients
print("\nEquation with numerical coefficients:")
print(f"F(theta) = {coeff_A0} + ({coeff_A1})*cos(theta) + ({coeff_B1})*sin(theta) + ({coeff_A2})*cos(2*theta) + ({coeff_B2})*sin(2*theta)")

# This illustrates how the overall poloidal structure is a sum of simple sinusoidal harmonics,
# which is the fundamental concept of the Fourier series expansion technique.