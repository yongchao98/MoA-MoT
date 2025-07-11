import math

# The height is given in the problem as 350 * g meters.
h_coefficient = 350

# The derived acceleration is a = g/7, so the fraction is 1/7.
# The time equation is t = sqrt(2*h/a)
# t = sqrt(2 * (350*g) / (g/7))
# The 'g' terms cancel out, so we can calculate the result from the coefficients.
t_squared = 2 * h_coefficient * 7
time = math.sqrt(t_squared)

print(f"The equation for the square of the time is t^2 = 2 * h / a")
print(f"Substituting the values: t^2 = 2 * ({h_coefficient}*g) / (g/7)")
print(f"The parameter 'g' cancels out: t^2 = 2 * {h_coefficient} * 7")
print(f"t^2 = {int(t_squared)}")
print(f"The time taken is the square root of {int(t_squared)}, which is {int(time)} seconds.")

# Final answer expression
# The final result is a numerical value.
final_answer = int(time)