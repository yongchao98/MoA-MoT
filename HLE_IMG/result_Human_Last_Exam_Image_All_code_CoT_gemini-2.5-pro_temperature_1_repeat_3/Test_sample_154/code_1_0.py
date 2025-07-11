# Step 1: Define the constants from the power system diagram and problem description.

# Local Power Generation (in MW)
p_plant_ga_gb_gc_units = 2
p_plant_ga_gb_gc_rating = 180
num_ga_gb_gc_plants = 3

p_plant_gd_ge_gf_units = 3
p_plant_gd_ge_gf_rating = 15
num_gd_ge_gf_plants = 3

# Transformer Capacities (in MVA)
sub_c_units = 4
sub_c_rating = 150
num_sub_c = 3

sub_d_units = 3
sub_d_rating = 63

sub_e1_units = 3
sub_e1_rating = 50

sub_e2_units = 3
sub_e2_rating = 40

# Loss Percentages
base_loss_rate = 0.02  # 2%
# From option D, the harmonic resonance increases losses by 8.5%
harmonic_loss_increase_rate = 0.085 # 8.5%


# Step 2: Calculate the total local power generation (P_local).
p_local_ga_gb_gc = num_ga_gb_gc_plants * p_plant_ga_gb_gc_units * p_plant_ga_gb_gc_rating
p_local_gd_ge_gf = num_gd_ge_gf_plants * p_plant_gd_ge_gf_units * p_plant_gd_ge_gf_rating
p_local_total = p_local_ga_gb_gc + p_local_gd_ge_gf

print("Calculation Steps:")
print(f"1. Total Local Generation (P_local):")
print(f"   - Power from GA, GB, GC plants = {num_ga_gb_gc_plants} * ({p_plant_ga_gb_gc_units} * {p_plant_ga_gb_gc_rating} MW) = {p_local_ga_gb_gc} MW")
print(f"   - Power from GD, GE, GF plants = {num_gd_ge_gf_plants} * ({p_plant_gd_ge_gf_units} * {p_plant_gd_ge_gf_rating} MW) = {p_local_gd_ge_gf} MW")
print(f"   - P_local = {p_local_ga_gb_gc} + {p_local_gd_ge_gf} = {p_local_total} MW\n")


# Step 3: Estimate the total system load (P_load) based on transformer capacities.
s_total_sub_c = num_sub_c * sub_c_units * sub_c_rating
s_total_sub_d = sub_d_units * sub_d_rating
s_total_sub_e1 = sub_e1_units * sub_e1_rating
s_total_sub_e2 = sub_e2_units * sub_e2_rating
s_total_transformers = s_total_sub_c + s_total_sub_d + s_total_sub_e1 + s_total_sub_e2

# Assumption: P_load (MW) is numerically equal to the total transformer capacity (MVA)
p_load_assumed = s_total_transformers

print(f"2. Estimated System Load (P_load):")
print(f"   - Total transformer capacity = ({num_sub_c}*{sub_c_units}*{sub_c_rating}) + ({sub_d_units}*{sub_d_rating}) + ({sub_e1_units}*{sub_e1_rating}) + ({sub_e2_units}*{sub_e2_rating}) = {s_total_transformers} MVA")
print(f"   - Assuming P_load is numerically equal to total MVA, P_load = {p_load_assumed} MW\n")


# Step 4: Calculate total system losses (P_losses).
total_loss_rate = base_loss_rate + harmonic_loss_increase_rate
p_losses = p_load_assumed * total_loss_rate

print(f"3. Total System Losses (P_losses):")
print(f"   - Base loss rate = {base_loss_rate*100}%")
print(f"   - Harmonic loss increase (from option D) = {harmonic_loss_increase_rate*100}%")
print(f"   - Total loss rate = {base_loss_rate*100}% + {harmonic_loss_increase_rate*100}% = {total_loss_rate*100:.1f}%")
print(f"   - P_losses = {total_loss_rate:.3f} * {p_load_assumed} MW = {p_losses:.2f} MW\n")


# Step 5: Calculate the power supplied by the external network (P_ext).
# P_ext = P_load + P_losses - P_local
p_ext_calculated = p_load_assumed + p_losses - p_local_total

print(f"4. Power from External Network (P_ext):")
print(f"   - P_ext = P_load + P_losses - P_local")
print(f"   - P_ext = {p_load_assumed} MW + {p_losses:.2f} MW - {p_local_total} MW = {p_ext_calculated:.2f} MW\n")

# Step 6: Compare with the answer choices.
print("5. Conclusion:")
print(f"The calculated external power supply is {p_ext_calculated:.2f} MW.")
print("This value is closest to 1273.2 MW from option D.")
print("The calculation used the 8.5% loss increase from option D, making it the most consistent choice.")
