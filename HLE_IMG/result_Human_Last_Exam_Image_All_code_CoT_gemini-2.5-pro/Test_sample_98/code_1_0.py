import math

# The Tarsus Medallion is widely dated to the reign of Antiochus IV Epiphanes,
# which began in 175 B.C. We will use this year to calculate the century.
year_bc = 175

# The formula to calculate the century for a B.C. year is: floor((year - 1) / 100) + 1
# Let's show the calculation with all the numbers.
term1 = year_bc
term2 = 1
term3 = 100
term4 = 1

century = math.floor((term1 - term2) / term3) + term4

print(f"To find the century for the year {year_bc} B.C., we use the formula: floor(({term1} - {term2}) / {term3}) + {term4}")
print(f"The calculation resolves to: floor(({term1 - term2}) / {term3}) + {term4}")
print(f"Which simplifies to: floor({(term1 - term2) / term3}) + {term4}")
print(f"This gives us: {math.floor((term1 - term2) / term3)} + {term4}")
print(f"Resulting in the century: {century}")

# Output the final answer in a readable format.
# For the number 2, the ordinal suffix is 'nd'.
print(f"\nTherefore, archaeologists ascribe this medallion to the {century}nd century B.C.")