# The medallion depicts the Sasanian King Bahram II.
# His reign started in the year 274 A.D.
year = 274

# To calculate the century from a year A.D.:
# 1. Subtract 1 from the year.
# 2. Perform integer division by 100.
# 3. Add 1 to the result.
# Equation: (year - 1) // 100 + 1
century = (year - 1) // 100 + 1
num1 = year
num2 = 1
num3 = 100
num4 = 1

# We need to output each number in the final equation.
# (274 - 1) // 100 + 1 = 3
print(f"Using the reign's start year, the calculation is:")
print(f"({num1} - {num2}) // {num3} + {num4} = {century}")
print("\nTherefore, the medallion is from the:")
print(f"{century}rd century A.D.")
