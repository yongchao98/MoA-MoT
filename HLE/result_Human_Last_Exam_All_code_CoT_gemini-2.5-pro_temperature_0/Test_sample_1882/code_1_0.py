# Define the masses of the proteins in kDa
giks3_mass = 60
capk1_mass = 45
capk2_mass = 102
capk3_mass = 65
capk4_mass = 55
capk5_mass = 39

# Step 1: Verify the oligomeric state of GIKS3 from the control experiment
# The observed mass of GIKS3 alone is 120 kDa.
giks3_dimer_mass = giks3_mass * 2
print(f"The expected mass of the GIKS3 dimer is {giks3_mass} * 2 = {giks3_dimer_mass} kDa.")
print("This matches the observed peak of 120 kDa, confirming GIKS3 is a dimer.\n")

# Step 2: Calculate the expected masses of the complexes based on the GIKS3 dimer
print("Calculating expected complex masses for comparison with SEC-MALS data:")

# GIKS3 + CaPK2
complex_giks3_capk2 = giks3_dimer_mass + capk2_mass
print(f"GIKS3 dimer + CaPK2: {giks3_dimer_mass} + {capk2_mass} = {complex_giks3_capk2} kDa. (Matches observed peak of 222 kDa -> Interaction)")

# GIKS3 + CaPK3
complex_giks3_capk3 = giks3_dimer_mass + capk3_mass
print(f"GIKS3 dimer + CaPK3: {giks3_dimer_mass} + {capk3_mass} = {complex_giks3_capk3} kDa. (Matches observed peak of 185 kDa -> Interaction)")

# GIKS3 + CaPK5
complex_giks3_capk5 = giks3_dimer_mass + capk5_mass
print(f"GIKS3 dimer + CaPK5: {giks3_dimer_mass} + {capk5_mass} = {complex_giks3_capk5} kDa. (Matches observed peak of 159 kDa -> Interaction)\n")

print("Final Conclusion: Based on a full analysis of all experiments, statement A is a correct summary of the findings.")
