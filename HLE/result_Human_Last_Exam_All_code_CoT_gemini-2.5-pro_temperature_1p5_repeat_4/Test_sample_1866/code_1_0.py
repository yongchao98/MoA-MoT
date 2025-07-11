# Step 1: Calculate cuts for the X-axis
# To cut a 4cm dimension requires a middle cut and quarter cuts.
# The middle cut is on a 4cm thick piece, requiring ceil(4/2) = 2 passes.
cuts_x_mid = 2
# The quarter cuts are on 2cm thick pieces, requiring ceil(2/2) = 1 pass.
cuts_x_quarter = 1
total_cuts_x = cuts_x_mid + cuts_x_quarter

# Step 2: Calculate cuts for the Z-axis
# After X-cuts, pieces have a 1cm dimension. We arrange them so the horizontal
# knife path is 2cm (by making a 2-piece-wide stack).
# The middle cut on the 4cm length requires ceil(2/2) = 1 pass.
cuts_z_mid = 1
# The quarter cuts on the remaining 2cm length require ceil(2/2) = 1 pass.
cuts_z_quarter = 1
total_cuts_z = cuts_z_mid + cuts_z_quarter

# Step 3: Calculate cuts for the Y-axis
# After X and Z cuts, pieces have a 1cm height.
# The middle cut on the 4cm length requires ceil(1/2) = 1 pass.
cuts_y_mid = 1
# The quarter cuts on the remaining 2cm length require ceil(1/2) = 1 pass.
cuts_y_quarter = 1
total_cuts_y = cuts_y_mid + cuts_y_quarter

# Step 4: Calculate the total minimum cuts
total_cuts = total_cuts_x + total_cuts_z + total_cuts_y

print("Calculation of the minimum number of cuts:")
print(f"Cuts for X-axis = {cuts_x_mid} (mid) + {cuts_x_quarter} (quarters) = {total_cuts_x}")
print(f"Cuts for Z-axis = {cuts_z_mid} (mid) + {cuts_z_quarter} (quarters) = {total_cuts_z}")
print(f"Cuts for Y-axis = {cuts_y_mid} (mid) + {cuts_y_quarter} (quarters) = {total_cuts_y}")
print("-" * 20)
print(f"Total Minimum Cuts = {total_cuts_x} + {total_cuts_z} + {total_cuts_y} = {total_cuts}")
