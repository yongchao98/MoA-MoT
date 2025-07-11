# The Mausoleum of Augustus was intended to be a lasting testament to the emperor's legacy.
# The most significant part of this legacy was described in the Res Gestae Divi Augusti,
# the "Accomplishments of the Divine Augustus," which was posted at the entrance.

# This script performs a calculation based on facts related to that document.

# 1. Key numbers are defined based on historical facts.
#    - The Res Gestae is composed of a preamble (1), 35 main sections, and an addendum (1).
#    - Section 15 of the Res Gestae mentions three occasions of cash donations to the plebs.
#    - The text was inscribed on two bronze pillars.
res_gestae_parts = 1 + 35 + 1
donation_occasions = 3
bronze_pillars = 2

# 2. Perform the calculation.
result = (res_gestae_parts * donation_occasions) - bronze_pillars

# 3. Print the equation and its components for clarity.
print(f"The number of total parts in the Res Gestae Divi Augusti is {res_gestae_parts}.")
print(f"The number of cash donation occasions noted in section 15 is {donation_occasions}.")
print(f"The number of bronze pillars the text was inscribed on is {bronze_pillars}.")
print("\nThe final calculation is:")
print(f"({res_gestae_parts} * {donation_occasions}) - {bronze_pillars} = {result}")
