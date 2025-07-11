import numpy as np

# Poles of B(z)
# From P1(z) = 4*z**4 - z**3 + z**2 + 1 = 0
# The sum of roots is -(-1)/4 = 1/4
sum_poles_B1 = 0.25
num_poles_B1 = 4

# From P2(z) = z**4 + z**2 - z + 4 = 0
# The sum of roots is -(0)/1 = 0
sum_poles_B2 = 0
num_poles_B2 = 4

# Total poles for B(z)
sum_poles_B = sum_poles_B1 + sum_poles_B2
num_poles_B = num_poles_B1 + num_poles_B2

# Poles of E(z) are hypothesized to be at z=1 and z=2
poles_E = [1, 2]
sum_poles_E = sum(poles_E)
num_poles_E = len(poles_E)

# Total poles for S(z) = E(z)B(z)
total_sum_of_poles = sum_poles_B + sum_poles_E
total_num_of_poles = num_poles_B + num_poles_E

# Calculate the average
average_z = total_sum_of_poles / total_num_of_poles

# Output the result
print(f"The sum of poles from the first polynomial for B(z) is: {sum_poles_B1}")
print(f"The number of poles from the first polynomial for B(z) is: {num_poles_B1}")
print(f"The sum of poles from the second polynomial for B(z) is: {sum_poles_B2}")
print(f"The number of poles from the second polynomial for B(z) is: {num_poles_B2}")
print(f"The sum of poles for E(z) is: {sum_poles_E}")
print(f"The number of poles for E(z) is: {num_poles_E}")
print(f"The total sum of all poles is: {total_sum_of_poles}")
print(f"The total number of poles is: {total_num_of_poles}")
print(f"The average value of the complex coordinates is: {total_sum_of_poles} / {total_num_of_poles} = {average_z}")
print(f"Final Answer: The average value is {average_z}")