# The medallion is generally dated by archaeologists to the year 215 A.D.
year = 215
one = 1
one_hundred = 100

# To find the century from the year, we use the formula: (year - 1) // 100 + 1
century = (year - one) // one_hundred + one

# We are asked to show the numbers in the final equation.
# The result will be a number, which we will format as an ordinal (e.g., 3rd).
ordinal_suffix = "rd" # for 3

# Print the final answer.
print(f"The equation to find the century from the year {year} is: ({year} - {one}) // {one_hundred} + {one} = {century}")
print(f"Therefore, archaeologists ascribe this medallion to the {century}{ordinal_suffix} century A.D.")