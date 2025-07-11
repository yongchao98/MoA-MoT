# Step 1: Define the atomic composition of the neutral 1,3,6,8-tetramethoxyacridine core.
# The acridine skeleton is C13H9N.
# Four methoxy groups (-OCH3) are added, replacing four H atoms.
c_acridine_skel = 13
h_acridine_skel = 9
n_acridine_skel = 1
o_acridine_skel = 0

num_methoxy_groups = 4
c_per_methoxy = 1
h_per_methoxy = 3
o_per_methoxy = 1

c_core = c_acridine_skel + num_methoxy_groups * c_per_methoxy
h_core = (h_acridine_skel - num_methoxy_groups) + (num_methoxy_groups * h_per_methoxy)
n_core = n_acridine_skel
o_core = o_acridine_skel + num_methoxy_groups * o_per_methoxy

# Step 2: Define the atomic composition of the N-substituent from methyl-3-aminopropionate.
# The substituent is -CH2-CH2-COOCH3.
c_substituent = 1 + 1 + 1 + 1
h_substituent = 2 + 2 + 3
n_substituent = 0
o_substituent = 2

# Step 3: Calculate the total atomic composition of the cation B.
c_total = c_core + c_substituent
h_total = h_core + h_substituent
n_total = n_core + n_substituent
o_total = o_core + o_substituent

# Step 4: Print the detailed calculation and the final molecular formula.
print("Calculation of the molecular formula for the cation of compound B:")
print("\n1. Atoms in the neutral 1,3,6,8-tetramethoxyacridine core:")
print(f"   C = {c_acridine_skel} + {num_methoxy_groups} * {c_per_methoxy} = {c_core}")
print(f"   H = ({h_acridine_skel} - {num_methoxy_groups}) + {num_methoxy_groups} * {h_per_methoxy} = {h_core}")
print(f"   N = {n_core}")
print(f"   O = {num_methoxy_groups} * {o_per_methoxy} = {o_core}")

print("\n2. Atoms in the N-substituent (-CH2CH2COOCH3):")
print(f"   C = {c_substituent}")
print(f"   H = {h_substituent}")
print(f"   O = {o_substituent}")

print("\n3. Total atoms in the final cation B+:")
print(f"   C = {c_core} + {c_substituent} = {c_total}")
print(f"   H = {h_core} + {h_substituent} = {h_total}")
print(f"   N = {n_total}")
print(f"   O = {o_core} + {o_substituent} = {o_total}")

print(f"\nThe molecular formula of compound B (cation part) is: C{c_total}H{h_total}NO{o_total}")