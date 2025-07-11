# Define the total number of chairs.
total_chairs = 20

# The optimal strategy is to seat people non-adjacently, creating a pattern
# like (Person, Empty, Person, Empty, ...).
# This means for every two chairs, one can be occupied.
# To find the maximum number of occupied chairs, we divide the total number of chairs by 2.
max_occupied = total_chairs // 2

# Print the final equation showing the calculation as requested.
print(f"To find the maximum number of occupied chairs, we follow a strategy where no two people sit next to each other.")
print(f"This means for {total_chairs} chairs, we can occupy every second chair.")
print(f"The final calculation is:")
print(f"{total_chairs} / 2 = {max_occupied}")