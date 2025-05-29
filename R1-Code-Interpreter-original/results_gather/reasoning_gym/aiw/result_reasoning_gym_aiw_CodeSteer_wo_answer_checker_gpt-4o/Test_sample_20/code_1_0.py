# Define the sets for female colleagues in each circle
female_colleagues_mary = {"F1", "F2", "F3"}  # 3 female colleagues in Mary's circle
female_colleagues_thomas = {"F4"}  # 1 female colleague in Thomas's circle

# Matilda is common in both circles
matilda = {"Matilda"}

# Calculate the total unique female colleagues Matilda has
# Union of both sets minus Matilda
unique_female_colleagues = (female_colleagues_mary | female_colleagues_thomas | matilda) - matilda

# Output the result
print(f"<<<{len(unique_female_colleagues)}>>>")