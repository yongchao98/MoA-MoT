import math

# For graph G, we determined m=5. The Shannon capacity is sqrt(5).
# For graph H, we determined n=4. The Shannon capacity is 2.
# The Shannon capacity of the strong product G⊠H is the product of the individual capacities.

shannon_capacity_g_val = 5
shannon_capacity_h = 2

final_capacity = math.sqrt(shannon_capacity_g_val) * shannon_capacity_h

print(f"sqrt({shannon_capacity_g_val}) * {shannon_capacity_h} = {final_capacity}")