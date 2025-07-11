# This script determines the IUPAC name of the product from the described reaction.

# The reaction is a stereoselective Ireland-Claisen rearrangement.
# The product is a secondary amide with a new stereocenter at the alpha-carbon.
# The IUPAC name is constructed based on the structure of the major diastereomer formed.
# The name includes all the locant numbers (2, 5, 1) as requested.
#
# (2S) - Absolute configuration of the new stereocenter at position 2.
# 2 - The locant for the main substituent on the propanamide chain.
# [((S)-5-methylcyclopent-1-en-1-yl)methyl] - The substituent group.
#     (S) - The existing stereocenter on the cyclopentene ring.
#     5 - Locant for the methyl group.
#     1 - Locant for the double bond ('en').
#     1 - Locant for the attachment point ('yl').
# N - Indicates the substituent is on the Nitrogen atom.
# ((S)-1-phenylethyl) - The substituent on the Nitrogen.
#     (S) - The stereocenter of the auxiliary.
#     1 - Locant for the phenyl group.
# propanamide - The parent amide chain.

final_product_iupac_name = "(2S)-2-[((S)-5-methylcyclopent-1-en-1-yl)methyl]-N-((S)-1-phenylethyl)propanamide"

print(final_product_iupac_name)