# Define the properties of the original nanocar chassis
chassis_formula = "C44H26O4"
num_hydrogens_in_chassis = 26

print("This program calculates the number of fluorine atoms in a hypothetical 'perfluoronanocar'.")
print("We base our calculation on the chassis of the first nanocar, synthesized by the Tour group (2005).")
print(f"The molecular formula of this chassis is {chassis_formula}.")
print("-" * 60)

print(f"The term 'perfluoro-' means every hydrogen atom (H) is replaced by a fluorine atom (F).")
print(f"The nanocar chassis has {num_hydrogens_in_chassis} hydrogen atoms available for replacement.")

# The calculation is a 1-to-1 replacement
fluorine_per_hydrogen = 1
total_fluorine_atoms = fluorine_per_hydrogen * num_hydrogens_in_chassis

print("\nThe final calculation is based on this one-to-one replacement:")
print(f"  {num_hydrogens_in_chassis} (hydrogen atoms in original molecule)")
print(f"x {fluorine_per_hydrogen} (fluorine atom per hydrogen atom)")
print("---------------------------------------")
print(f"= {total_fluorine_atoms} (total fluorine atoms in perfluoronanocar)")
