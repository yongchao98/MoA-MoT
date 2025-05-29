# Define the social circles using sets
sarah_circle = {"Sarah", "Male1", "Male2", "Male3", "Female1"}
john_circle = {"John", "Male4", "Male5", "Male6", "Female2", "Matilda"}

# Identify Matilda as the common colleague
common_colleagues = sarah_circle.intersection(john_circle)

# Since we know Matilda is the only common colleague, we can remove her from the set
common_colleagues.discard("Matilda")

# Calculate the unique female colleagues Matilda has
# Exclude Matilda herself from the count
sarah_females = {"Female1"}
john_females = {"Female2", "Matilda"}

# Use set union to find all unique female colleagues
unique_female_colleagues = sarah_females.union(john_females) - {"Matilda"}

# Output the final answer
print(f"<<<{len(unique_female_colleagues)}>>>")