# Evaluate each statement about North American Cordilleran geology
# C = Consensus, D = Debated

# (1) The Morrison formation represents a foredeep deposit.
# Status: D (Debated). The term "foredeep" is a specific classification that is debated for the Morrison Formation.
s1 = "D"

# (2) Metamorphic core complexes formed in response to a slab window.
# Status: D (Debated). A slab window is one of several competing hypotheses for the formation of metamorphic core complexes.
s2 = "D"

# (3) An orogenic plateau (the "Nevadaplano") existed in present-day Nevada by the end of the Cretaceous.
# Status: D (Debated). The existence and elevation of the Nevadaplano is a major area of active debate.
s3 = "D"

# (4) An orogenic plateau in the Arizona portion of the Cordillera (the "Arizonaplano") was constructed via magmatic thickening of the crust in the late Cretaceous.
# Status: D (Debated). Similar to the Nevadaplano, the Arizonaplano hypothesis is highly debated.
s4 = "D"

# (5) Laramide-style structures are generally inboard of Sevier-style structures.
# Status: C (Consensus). This is a fundamental and widely accepted spatial relationship in North American tectonics.
s5 = "C"

# (6) The Sierra Nevada batholith and Idaho Batholith were formed by the ancestral Cordilleran arc.
# Status: C (Consensus). This is the textbook, consensus model for the origin of these batholiths.
s6 = "C"

# (7) A series of large-magnitude ignimbrite eruptions began in Idaho and Montana in the Eocene and propagated generally to the southwest.
# Status: C (Consensus). This accurately describes the well-established track of the Yellowstone hotspot.
s7 = "C"

# (8) The Colorado Plateau reached modern elevations by the late Cretaceous, and may have been uplifted further by Cenozoic delamination of mantle lithosphere.
# Status: D (Debated). The timing and mechanism of the Colorado Plateau's uplift are classic, unresolved debates in geology.
s8 = "D"

# (9) The Cascades arc reached modern elevation by the late Oligocene, but experienced a Pliocene-present flare-up.
# Status: D (Debated). The modern high peaks of the Cascades are much younger than the Oligocene, making the first part of the statement incorrect/non-consensus.
s9 = "D"

# (10) Laramide-style deformation was caused by subduction of the Shatsky conjugate.
# Status: D (Debated). This is one of several prominent, competing hypotheses for the cause of the Laramide orogeny.
s10 = "D"

# Combine the results into a single string
final_string = s1 + s2 + s3 + s4 + s5 + s6 + s7 + s8 + s9 + s10

print(final_string)