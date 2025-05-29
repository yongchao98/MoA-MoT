# Initial cost of the plane
plane_cost = 150000

# Monthly costs
hanger_rent = 8578246
fuel_cost = 2 * hanger_rent

# Total costs for 12 months
total_hanger_rent = hanger_rent * 12
total_fuel_cost = fuel_cost * 12

# Total cost for the first year
total_cost = plane_cost + total_hanger_rent + total_fuel_cost

print(total_cost)