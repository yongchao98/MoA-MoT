# The mapping class is (Da o Db)^9. We are interested in the total power.
total_power = 9

# For a torus with one boundary component, there is a fundamental relation:
# (Da o Db)^3 = Dd
# where Dd is a Dehn twist about the boundary.
# The power in this relation is 3.
relation_power = 3

# The fractional Dehn twist coefficient can be found by expressing the
# original mapping class in terms of Dd.
# (Da o Db)^9 = ((Da o Db)^3)^3 = (Dd)^3
# The coefficient is the final exponent on Dd.
# This can be calculated by dividing the total power by the relation power.
coefficient = total_power / relation_power

# Print the numbers used in the final calculation and the result.
print(f"{total_power} / {relation_power} = {int(coefficient)}")