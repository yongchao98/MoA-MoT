# The plan to derive the magnetization curve is as follows:
# 1. Model the superconductor as a thin slab of width 2a.
# 2. Apply the Bean critical-state model, where induced current density is +/- Jc.
# 3. Relate the field penetration depth 'd' to the applied field 'H' using Ampere's law, yielding d = H / Jc.
# 4. Calculate the internal magnetic field B_ind(x) produced by these screening currents.
# 5. Average B_ind(x) over the slab's cross-section to find the magnetization M.
# 6. Substitute d = H / Jc to get the final expression M(H).

print("This script presents the derived analytical expression for the initial magnetization of a superconducting slab.")
print("\n" + "="*70 + "\n")

# Define the symbols used in the equation for clarity
print("Symbol Definitions:")
print("  M: Magnetization")
print("  H: Applied magnetic field")
print("  Jc: Critical current density")
print("  a: Half-width of the superconducting slab")
print("-" * 70)

# Present the final derived formula
# This formula is valid for 0 <= H <= Hp, where the full penetration field Hp = Jc * a.
final_formula_v1 = "M(H) = -H * (1 - H / (2 * Jc * a))"
final_formula_v2 = "M(H) = -H + H**2 / (2 * Jc * a)"
print("The analytical expression for the initial magnetization curve is:")
print(f"  {final_formula_v1}")
print("\nAlternatively, in expanded form:")
print(f"  {final_formula_v2}\n")

# Fulfill the requirement to "output each number in the final equation"
# by breaking it down into a canonical form and listing the coefficients.
print("To highlight each number in the final equation, let's examine the expanded form:")
print(f"  M(H) = (-1) * H**1 + (1 * H**2) / (2 * Jc * a)\n")

print("The constant numbers (coefficients and powers) in this expression are:")
print(f"  Coefficient of the first term: -1")
print(f"  Power of H in the first term:   1")
print(f"  Coefficient of the second term: 1")
print(f"  Power of H in the second term:  2")
print(f"  Constant in the denominator:    2")