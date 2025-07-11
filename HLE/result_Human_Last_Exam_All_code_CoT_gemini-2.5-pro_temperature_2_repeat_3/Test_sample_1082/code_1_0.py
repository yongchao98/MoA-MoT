# All values are in dollars.

# Revenue from selling the final product.
revenue = 50

# Environmental cost for water used in the process.
water_cost = 10

# Environmental cost for energy used in the process.
energy_cost = 15

# Calculate the total environmental costs.
total_environmental_cost = water_cost + energy_cost

# Calculate the Sustainable Value Added (SVA).
sva = revenue - total_environmental_cost

# The problem requires printing the final equation with each number.
# The equation for SVA is Revenue - (Water Cost + Energy Cost).
print(f"Sustainable Value Added (SVA) Calculation:")
print(f"SVA = {revenue} - ({water_cost} + {energy_cost})")
print(f"SVA = {revenue} - {total_environmental_cost}")
print(f"Final SVA = ${sva}")