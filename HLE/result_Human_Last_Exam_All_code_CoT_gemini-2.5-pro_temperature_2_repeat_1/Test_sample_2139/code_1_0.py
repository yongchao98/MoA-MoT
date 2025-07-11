import numpy as np
from scipy.special import gamma

# Calculate the constant C which is y(0)
# y(0) = (128 * 3^(1/6) * Gamma(2/3))^(-1)
val_128 = 128
val_cbrt3 = np.power(3, 1/6)
val_gamma = gamma(2/3)
y0 = (val_128 * val_cbrt3 * val_gamma)**-1

print(f"The equation for the radius at t=pi/4 is:")
print(f"y(pi/4) = y(0) * ( (pi/4)^4 + 1 )^(-1/2) * ( cos(pi/4) )^(1/4)")
print(f"Where:")
print(f"y(0) = (128 * 3^(1/6) * Gamma(2/3))^(-1) = ({val_128} * {val_cbrt3:.6f} * {val_gamma:.6f})^-1 = {y0:.6f}")

# Calculate the second term at t = pi/4
t = np.pi/4
term2_base = (t**4 + 1)
term2 = term2_base**(-0.5)
print(f"( (pi/4)^4 + 1 )^(-1/2) = ( ({t:.6f})^4 + 1 )^(-1/2) = {term2:.6f}")


# Calculate the third term at t = pi/4
term3_base = np.cos(t)
term3 = term3_base**(0.25)
print(f"( cos(pi/4) )^(1/4) = ( {term3_base:.6f} )^(1/4) = {term3:.6f}")

# Calculate the final result
y_pi_over_4 = y0 * term2 * term3

print(f"\nFinal calculation:")
print(f"y(pi/4) = {y0:.6f} * {term2:.6f} * {term3:.6f}")
print(f"y(pi/4) = {y_pi_over_4}")
