# Step-by-step analysis to determine the IUPAC name of C7H14

# 1. Analyze the Molecular Formula
# Formula: C7H14
# The general formula for an acyclic alkane is C(n)H(2n+2). For n=7, this would be C7H16.
# The given formula C7H14 has two fewer hydrogens than the corresponding saturated alkane.
# This indicates one degree of unsaturation.
molecular_formula = "C7H14"
carbon_atoms = 7
hydrogen_atoms = 14
degrees_of_unsaturation = (2 * carbon_atoms + 2 - hydrogen_atoms) / 2

print("Step 1: Analysis of the Molecular Formula")
print(f"The molecular formula is {molecular_formula}.")
print(f"This corresponds to a degree of unsaturation of {int(degrees_of_unsaturation)}, indicating one C=C double bond or one ring.")
print("-" * 30)

# 2. Analyze the 13C NMR Data
# Data: 145(s), 112(t), 48(t), 27(d), 22(q), 21(q).
# Note on multiplicity: (s)=quaternary C, (d)=CH, (t)=CH2, (q)=CH3.
# The chemical shifts at 145 ppm and 112 ppm are in the alkene region (100-150 ppm). This confirms the C=C double bond.
# - 145(s): A quaternary (s) sp2 carbon (no attached H), part of the double bond.
# - 112(t): A secondary (t) sp2 carbon (a CH2 group), part of the double bond.
# These two signals confirm a terminal, disubstituted alkene group: >C=CH2.

print("Step 2: Analysis of 13C NMR Data")
print("The chemical shifts at 145 ppm and 112 ppm confirm a C=C double bond.")
print("The multiplicities reveal the following alkene structure fragment:")
print(" - A carbon with signal 145(s) is a quaternary C in the double bond.")
print(" - A carbon with signal 112(t) is a CH2 group in the double bond.")
print("This gives the core structure: >C=CH2")
print("-" * 30)


# 3. Assemble the Structure Fragments
# The molecule has 7 carbons, but only 6 NMR signals are given, confirming signal overlap.
# Let's account for all atoms from the signals:
# C from 145(s): 1C, 0H
# CH2 from 112(t): 1C, 2H
# CH2 from 48(t): 1C, 2H
# CH from 27(d): 1C, 1H
# CH3 from 22(q): 1C, 3H
# CH3 from 21(q): 1C, 3H
# Total accounted for from 6 signals: 6 Carbons, 11 Hydrogens.
# The missing atoms from C7H14 are 1 Carbon and 3 Hydrogens, which is a CH3 group.
# This means one of the methyl (q) signals represents two equivalent CH3 groups.
# The combination of a CH signal at 27(d) and two equivalent CH3 groups (from the 21q signal) strongly suggests an isopropyl group: -CH(CH3)2.
# The remaining alkyl fragments are a CH2 group at 48(t) and a CH3 group at 22(q). These form an ethyl group: -CH2-CH3.
# The full structure is formed by attaching both the ethyl group and the isopropyl group to the quaternary carbon of the >C=CH2 core.

print("Step 3: Assembling the Fragments")
print("The alkyl signals are: 48(t), 27(d), 22(q), 21(q).")
print("A CH group at 27(d) and an overlapping signal for two CH3 groups at 21(q) indicate an isopropyl group: -CH(CH3)2.")
print("A CH2 group at 48(t) and a CH3 group at 22(q) indicate an ethyl group: -CH2CH3.")
print("Attaching both groups to the quaternary carbon of the >C=CH2 core gives the following structure:")
print("      CH(CH3)2")
print("     /")
print("CH2=C")
print("     \\")
print("      CH2CH3")
print("-" * 30)

# 4. Determine the IUPAC Name
# Structure: CH2=C(CH2CH3)(CH(CH3)2)
# To find the IUPAC name, we identify the longest carbon chain that includes the double bond and has the most substituents.
# - Chain via isopropyl group: Gives a butene parent with two substituents (ethyl at C2, methyl at C3).
# - Chain via ethyl group: Gives a butene parent with one substituent (isopropyl at C2).
# We choose the chain with more substituents.
# Numbering from the CH2= end gives the double bond the #1 position.
#       CH3 (substituent on C3)
#      |
# CH3--CH---C(CH2CH3)==CH2
# 4    3     2           1

print("Step 4: Determine the IUPAC Name")
print("The structure is identified as CH2=C(CH2CH3)(CH(CH3)2).")
print("The parent chain is a butene chain, numbered to give the double bond position 1 (1-butene).")
print("There is an ethyl group substituent at carbon 2.")
print("There is a methyl group substituent at carbon 3.")
final_name = "2-ethyl-3-methyl-1-butene"
print(f"\nThe final IUPAC name is: {final_name}")