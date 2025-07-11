# Liquid Crystal Molecular Design Script

# 1. Define the design parameters based on the provided guidelines.
#    - The general structure is C(n)H(2n+1)-Ph-CN.
#    - 'Ph' represents the single benzene ring core.
#    - '-CN' is the polar group.
#    - We will start with an alkyl chain length of n=5, a common choice for
#      achieving room-temperature liquid crystal phases.

n = 5

# 2. Define the components of the molecule
core = "Ph"  # Represents a para-substituted benzene ring (C6H4)
polar_group = "CN"
alkyl_chain_prefix = "C(n)H(2*n+1)"

# 3. Calculate the specific alkyl chain formula based on n=5
#    This fulfills the requirement to output each number in the final equation.
print("Designing the molecule using the general structure: CnH2n+1-Ph-CN")
print("----------------------------------------------------------")
print(f"Step 1: Set the alkyl chain length 'n'.")
print(f"Let n = {n}")
print("")

print("Step 2: Calculate the alkyl chain formula using the equation: C(n)H(2*n + 1)")
num_hydrogens = 2 * n + 1
print(f"The calculation is: 2 * {n} + 1 = {num_hydrogens}")
alkyl_chain_formula = f"C{n}H{num_hydrogens}"
print(f"The resulting flexible alkyl chain is: {alkyl_chain_formula}")
print("")

# 4. Assemble and describe the final molecule
molecule_name = "4-pentylbenzonitrile"
final_structure = f"{alkyl_chain_formula}-{core}-{polar_group}"

print("Step 3: Assemble the final molecular structure.")
print(f"Final Structure: {final_structure}")
print(f"Common Name: {molecule_name}")
print("")

print("Conclusion:")
print(f"The molecule '{molecule_name}' is an excellent candidate because:")
print("- It is based on a single benzene ring.")
print("- The combination of a rigid core (Ph), a flexible alkyl chain, and a polar group (CN) promotes nematic/smectic phases.")
print("- With n=5, its nematic phase transition temperature is ~35°C, which is close to room temperature.")