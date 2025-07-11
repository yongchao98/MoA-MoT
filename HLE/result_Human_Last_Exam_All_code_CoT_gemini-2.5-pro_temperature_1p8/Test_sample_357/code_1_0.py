# Step 1: Define the parameters based on the problem description.
# In a competitive market for handyman labor, homeowners (like Chelsea) will
# have their profits competed down to zero. The homeowner's profit is v - p.
# So, for any contract (v, p) in equilibrium, we must have p = v.

# Step 2: Determine the contract for the low type (L).
# The low type's utility is u_L = p - v^2. With p = v, this becomes u_L = v - v^2.
# To maximize their utility, the low type chooses v to maximize this function.
# The maximum of v - v^2 occurs at v = 0.5.
v_L = 0.5
p_L = v_L  # From the zero-profit condition

# Check the low type's productivity constraint: v_L <= 20 - (1/3)*v_L
is_vL_feasible = v_L <= 20 - (1/3) * v_L

# Step 3: Determine the contract for the high type (H).
# The high type's utility is u_H = p - v. With p = v, this becomes u_H = v - v = 0.
# The high type is indifferent between any contract (v, v).
# For the equilibrium to be separating, the low type must prefer (v_L, p_L) over (v_H, p_H).
# u_L at (v_L, p_L) is 0.5 - 0.5^2 = 0.25.
# u_L at (v_H, p_H) would be v_H - v_H^2. The condition 0.25 >= v_H - v_H^2 holds for all v_H.
# To find a unique value for v_H, we apply the "no distortion at the top" principle.
# This means the high type should perform at their maximum capacity, limited by their productivity theta_H.
theta_H = 30
v_H = theta_H
p_H = v_H  # From the zero-profit condition

# Check the high type's productivity constraint: v_H <= theta_H
is_vH_feasible = v_H <= theta_H

# Step 4: Finalize and print the contracts.
# The separating equilibrium consists of two contracts: one for the low type and one for the high type.
contract_L = (v_L, p_L)
contract_H = (v_H, p_H)

print("The separating equilibrium involves the following pair of contracts ((v_L, p_L), (v_H, p_H)):")
# We explicitly print each number in the final equation as requested.
print(f"(({contract_L[0]}, {contract_L[1]}), ({contract_H[0]}, {contract_H[1]}))")
