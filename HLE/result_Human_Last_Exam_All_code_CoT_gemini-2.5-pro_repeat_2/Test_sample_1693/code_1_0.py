# The task is to identify the largest French ship by displacement sunk by a U-boat
# before the armistice of June 22, 1940.

# Ship details
ship_name = "Champlain"
ship_type = "Passenger Liner"
# Gross Register Tonnage (GRT) is a common measure for civilian ships.
displacement_grt = 28124
sinking_date = "June 17, 1940"
cause = "Struck a mine laid by the German submarine U-65"

# Print the information clearly.
print(f"The largest French ship sunk by a U-boat before the 1940 armistice was the {ship_type}, '{ship_name}'.")
print(f"It was sunk on {sinking_date} after hitting a mine laid by U-65.")
print("The ship's size by Gross Register Tonnage (GRT) was:")
# The final part of the request is to output the number from the final "equation" or result.
print(displacement_grt)