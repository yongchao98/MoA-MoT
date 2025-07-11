# Plan:
# 1. Identify the rarest noble gas. Radon (Rn) is the rarest noble gas on Earth.
# 2. Define its abundance. We will use the mean atmospheric volume fraction, which is 6e-20.
#    This value is stable over human timescales, so the value for 2002 is the same.
# 3. Calculate the percentage by multiplying the fraction by 100.
# 4. Print the name of the gas, the final calculation showing all numbers, and the result.

# The rarest noble gas
rarest_gas_name = "Radon (Rn)"

# Its abundance as a volume fraction in the atmosphere
volume_fraction = 6e-20

# The number to multiply by to get a percentage
multiplier = 100

# The final calculation
percentage = volume_fraction * multiplier

print(f"The rarest noble gas on Earth is {rarest_gas_name}.")
print("Its percentage, based on its mean volume fraction in the atmosphere, is calculated below.")
print("\n--- Final Calculation ---")
print(f"Mean Volume Fraction: {volume_fraction}")
print(f"Multiplier: {multiplier}")
print(f"Final Equation: {volume_fraction} * {multiplier} = {percentage}")
print(f"-------------------------")
print(f"\nThe percentage of Radon on Earth is approximately: {percentage}%")