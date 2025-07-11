# This script determines the century of the medallion from the image.
# The medallion is identified as a piece from the Aboukir hoard,
# depicting Emperor Caracalla, who reigned from 198 to 217 A.D.
# We can use a year from his reign to find the century.

# Let's use a representative year from Caracalla's reign.
year = 215
period = "A.D."

# The formula to calculate the century from a year in the A.D. era is:
# century = (year // 100) + 1
# We will show the numbers used in this equation.

part1 = year // 100
part2 = 1
century = part1 + part2

print(f"The medallion dates to the reign of Emperor Caracalla (e.g., the year {year} {period}).")
print(f"To find the century, we perform the calculation: ({year} // 100) + {part2}")
print(f"The equation with its numbers is: {part1} + {part2} = {century}")
print(f"This means the medallion is from the 3rd century {period}.")
print(f"In the requested format, the answer is:")
print(f"{century} {period}")