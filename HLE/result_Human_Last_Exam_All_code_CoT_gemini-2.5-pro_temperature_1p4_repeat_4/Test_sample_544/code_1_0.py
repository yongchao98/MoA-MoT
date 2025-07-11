# Step 1: Define the reaction stoichiometry.
# The reaction consumes one equivalent of each of the three reactants to form
# one equivalent of the main product. These are the numbers in the balanced equation.
stoichiometry = {
    "Methyl phenyl sulfoxide": 1,
    "Triflic anhydride": 1,
    "Trimethylsilyl cyanide": 1,
    "Final Product": 1
}

print("Analysis of the reaction stoichiometry:")
print("The numbers from the balanced chemical equation are the stoichiometric coefficients:")
for component, number in stoichiometry.items():
    print(f"  {number} : {component}")


# Step 2: Determine the IUPAC name based on the product structure, C6H5-S-CH2-CN.
# The script will assemble the name from its constituent parts.

# The principal functional group is the nitrile (-CN). The parent chain containing this
# group is based on acetonitrile (CH3-CN).
parent_name = "acetonitrile"

# A substituent group, C6H5-S-, is attached to the acetonitrile skeleton.
# According to IUPAC nomenclature, this is named 'phenylsulfanyl'.
substituent_name = "phenylsulfanyl"

# The nitrile carbon is defined as carbon 1 (C1). Therefore, the carbon
# to which the substituent is attached is carbon 2 (C2). This is the number
# that will appear in the final name.
position = 2

# Step 3: Assemble the full IUPAC name.
# The format is position-(substituent_name)parent_name.
# Parentheses enclose the complex substituent name.
final_name = f"{position}-({substituent_name}){parent_name}"


# Step 4: Print the derived name and its components.
print("\nDerivation of the IUPAC name for the product (C6H5-S-CH2-CN):")
print(f"  - Parent Name (from -CH2-CN): {parent_name}")
print(f"  - Substituent Name (from C6H5-S-): {substituent_name}")
print(f"  - Position Number of Substituent: {position}")

print(f"\nThe final IUPAC name of the product is: {final_name}")