import numpy as np

N = 17
lamb = 2
a = 5
pi = np.pi

# The final expression is (N + lamb) * (y3_x0)**(lamb / a)
# y3_x0 was derived to be pi**4.5 / 320
y3_x0 = (pi**4.5) / 320

result = (N + lamb) * (y3_x0)**(lamb / a)

# To display the full equation as requested
# We compute each numerical part
part1 = N + lamb
part2 = y3_x0
part3 = lamb/a
print(f"({N} + {lamb}) * ({part2:.8f})^({lamb}/{a}) = {part1} * {part2**part3:.8f} = {result:.8f}")

# The final result in the requested format
# Since the problem expects a numerical answer as output.
print(f"The equation is: ({part1}) * ({part2})**({part3})")
print(f"The final calculated value is: {result}")
# I need to print the numbers in the final equation.
print(f"({(N + lamb)}) * (({pi}**4.5) / 320)**({lamb/a}) = {result}")