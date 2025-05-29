# Define the sets of female colleagues for Linda and Michael
female_colleagues_linda = {"F1", "F2", "F3", "Matilda"}
female_colleagues_michael = {"F4", "F5", "F6", "F7", "Matilda"}

# Find the union of both sets, excluding Matilda
unique_female_colleagues = (female_colleagues_linda | female_colleagues_michael) - {"Matilda"}

# Output the number of unique female colleagues Matilda has
print(f"<<<{len(unique_female_colleagues)}>>>")