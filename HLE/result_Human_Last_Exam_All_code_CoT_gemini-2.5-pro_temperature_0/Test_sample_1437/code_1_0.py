# Define the result from the theoretical derivation.
# In the static (w=0) and long-wavelength (q->0) limit, the ratio of the
# Lindhard function Pi(0,0) to the density of states at the Fermi level D(eps_F)
# is a constant for a 3D electron gas at T=0.
# The derived relationship is: Pi(0,0) = -D(eps_F).
result = -1

# --- Explanation and Calculation ---

print("The Lindhard polarization function, Pi(q, w), describes the charge density response of an electron gas to an external potential.")
print("We need to evaluate it in the static (w=0) and long-wavelength (momentum transfer q->0) limit.")
print("(Note: The prompt uses 'k' for momentum transfer, but we use the standard notation 'q' to avoid confusion with the electron wavevector k.)")

print("\nIn this limit, linear response theory shows that the function is equal to the negative of the density of states at the Fermi level, D(eps_F):")
print("Pi(q->0, w=0) = -D(eps_F)")

print("\nThe density of states D(eps_F) depends on material-specific properties like electron density. To find a universal numerical value, we evaluate the dimensionless Lindhard function, which is normalized by D(eps_F) itself.")

print("\nThe final calculation is as follows:")
print("Normalized Pi(0, 0) = Pi(0, 0) / D(eps_F)")
print("                   = -D(eps_F) / D(eps_F)")
print(f"                   = {result}")