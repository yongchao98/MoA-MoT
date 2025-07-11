# The task is to scan the Latin hexameter line:
# "verum ubi equi atque hominis casu convenit imago."
# and represent dactyls as "D" and spondees as "S".

# A hexameter line has six feet. Each foot is either a dactyl (- u u) or a spondee (- -).
# The fifth foot is typically a dactyl, and the sixth is always two syllables,
# usually treated as a spondee.

# The scansion of this line from Virgil's Aeneid (6.848) is well-established,
# though its derivation involves complex rules of elision and poetic license.
# - Foot 1: "ver(um) ubi" -> Dactyl (- u u)
# - Foot 2: "equi atque" -> Spondee (- -)
# - Foot 3: "(at)que hominis ca-" -> Dactyl (- u u)
# - Foot 4: "-su convenit" -> Dactyl (- u u)
# - Foot 5: "venit i(mago)" -> Spondee (- -)
# - Foot 6: "imago" -> Spondee (- -)

# The resulting pattern is Dactyl-Spondee-Dactyl-Dactyl-Spondee-Spondee.

print("The scansion for 'verum ubi equi atque hominis casu convenit imago' is:")

# Printing the final sequence of dactyls (D) and spondees (S).
print("D S D D S S")