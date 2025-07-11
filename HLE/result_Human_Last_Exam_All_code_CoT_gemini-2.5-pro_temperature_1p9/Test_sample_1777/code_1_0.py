# Define the given parameters
Eg = 3.0  # Band gap in eV
E1s = 1.0  # Energy of the 1s exciton peak in eV
n = 3     # Target principal quantum number

# The energy of a 2D exciton state is E_n = Eg - Ry* / (n - 0.5)^2
# First, calculate the effective Rydberg energy (Ry*) from the 1s state data.
# E1s = Eg - Ry* / (1 - 0.5)^2  =>  Ry* = (Eg - E1s) * (1 - 0.5)^2
Ry_star = (Eg - E1s) * (1 - 0.5)**2

# Now, calculate the energy for the n=3 state
# E3 = Eg - Ry* / (n - 0.5)^2
E3 = Eg - Ry_star / (n - 0.5)**2

# Print the final equation with all the numbers, as requested.
# The format is E(n) = Eg - Ry* / (n - 0.5)^2 = result
print(f"The energy for the n={n} state is calculated using the equation E(n) = Eg - Ry*/(n-0.5)^2:")
print(f"E({n}) = {Eg} - {Ry_star} / ({n} - 0.5)^2 = {E3:.2f} eV")