# Define the monoisotopic atomic masses
C_mass = 12.000000
H_mass = 1.007825
O_mass = 15.994915
N_mass = 14.003074

# Step 1: Calculate the mass of the itaconate adduct (C5H6O4) left after cleavage.
# The probe is derived from itaconic acid (C5H6O4) and propargylamine.
# After Michael addition and cleavage of the amide bond, the itaconate moiety remains.
num_C_adduct = 5
num_H_adduct = 6
num_O_adduct = 4

mass_itaconate_adduct = (num_C_adduct * C_mass) + (num_H_adduct * H_mass) + (num_O_adduct * O_mass)

# Step 2: Account for the additional mass modification.
# Formic acid can reduce the alkyne group on the probe before cleavage.
# The reduction of an alkyne (C≡C) to an alkane (CH2-CH2) adds 4 hydrogen atoms.
num_H_reduction = 4
mass_reduction = num_H_reduction * H_mass

# Step 3: Calculate the total variable modification mass (x).
# This is the sum of the itaconate adduct mass and the reduction mass.
total_modification_mass = mass_itaconate_adduct + mass_reduction

# Step 4: Print the final equation and the result
print("The mass calculation is based on the itaconate adduct plus a reduction.")
print(f"Mass of Itaconate Adduct (C5H6O4) = {mass_itaconate_adduct:.4f} Da")
print(f"Mass from reduction (4H) = {mass_reduction:.4f} Da")
print(f"Final Equation: x = Mass(C5H6O4) + Mass(4H)")
print(f"x = {mass_itaconate_adduct:.4f} + {mass_reduction:.4f}")
print(f"Calculated mass of modification (x) = {total_modification_mass:.4f} Da")
print(f"The closest integer value is {round(total_modification_mass)} Da.")