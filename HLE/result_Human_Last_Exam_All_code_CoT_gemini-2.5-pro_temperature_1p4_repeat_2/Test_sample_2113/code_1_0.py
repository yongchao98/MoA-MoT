import math

# Step 3: Compute the required values
# Calculate z1(pi/8)^2
# z1_at_pi_over_8_sq = cos^2(pi/8) = (1 + cos(pi/4))/2
cos_pi_over_4 = math.sqrt(2) / 2
z1_at_pi_over_8_sq = (1 + cos_pi_over_4) / 2

# Calculate y(pi/8)^2
# y(pi/8) = 2^(1/4), so y(pi/8)^2 = sqrt(2)
y_at_pi_over_8_sq = math.sqrt(2)

# Calculate the final ratio
final_value = z1_at_pi_over_8_sq / y_at_pi_over_8_sq

# We can also express the values symbolically for clarity in the output
# z1(pi/8)^2 = (sqrt(2) + 1) / (2 * sqrt(2))
z1_sq_symbolic_num = math.sqrt(2) + 1
z1_sq_symbolic_den = 2 * math.sqrt(2)

# y(pi/8)^2 = sqrt(2)
y_sq_symbolic = math.sqrt(2)

# Final result = z1_sq / y_sq = ( (sqrt(2)+1)/(2*sqrt(2)) ) / sqrt(2) = (sqrt(2)+1)/4
final_value_symbolic_num = math.sqrt(2) + 1
final_value_symbolic_den = 4

print(f"The trajectory z1(t) is cos(t). At t = pi/8, z1(pi/8) = cos(pi/8).")
print(f"The entanglement echo y(t) is sin(2t) / (cos(2t))^(3/2). At t = pi/8, y(pi/8) = 2^(1/4).")
print(f"We need to calculate (z1(pi/8) / y(pi/8))^2.")
print(f"The numerator is z1(pi/8)^2 = cos^2(pi/8) = (1 + cos(pi/4))/2 = ({z1_sq_symbolic_num / z1_sq_symbolic_den:.5f}) / ({y_sq_symbolic:.5f}) = {z1_at_pi_over_8_sq:.5f}")
print(f"The denominator is y(pi/8)^2 = (2^(1/4))^2 = sqrt(2) = {y_at_pi_over_8_sq:.5f}")
print(f"The final value is (z1(pi/8)^2) / (y(pi/8)^2) = {z1_at_pi_over_8_sq:.5f} / {y_at_pi_over_8_sq:.5f} = {final_value:.5f}")
print(f"The exact expression for the final value is (sqrt(2) + 1) / 4.")
