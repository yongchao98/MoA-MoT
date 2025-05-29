def calculate_female_colleagues_in_johns_circle(total_female_colleagues):
    """
    Calculate the number of female colleagues Matilda has in John's circle.
    
    :param total_female_colleagues: Total number of female colleagues in John's circle.
    :return: Number of female colleagues Matilda has, excluding herself.
    """
    # Matilda is one of the female colleagues, so subtract 1 to find her female colleagues
    return total_female_colleagues - 1

# Total female colleagues in John's circle
total_female_colleagues_john = 5

# Calculate the number of female colleagues Matilda has
female_colleagues_matilda = calculate_female_colleagues_in_johns_circle(total_female_colleagues_john)

# Output the result
print(f"<<<{female_colleagues_matilda}>>>")