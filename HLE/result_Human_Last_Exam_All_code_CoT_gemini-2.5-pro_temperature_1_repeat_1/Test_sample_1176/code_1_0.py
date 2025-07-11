# This script evaluates 10 statements about the geology of the North American Cordillera
# and determines whether each is a consensus (C) or debated (D) view.

# Statement (1): The Morrison formation represents a foredeep deposit.
# Verdict: C (Consensus). This is a textbook interpretation of a retroarc foreland basin.
s1 = "C"

# Statement (2): Metamorphic core complexes formed in response to a slab window.
# Verdict: D (Debated). The cause of extension is highly debated; a slab window is one of several hypotheses.
s2 = "D"

# Statement (3): An orogenic plateau ("Nevadaplano") existed in Nevada by the end of the Cretaceous.
# Verdict: D (Debated). The existence, scale, and timing of this proposed plateau are not settled.
s3 = "D"

# Statement (4): An orogenic plateau ("Arizonaplano") was constructed via magmatic thickening.
# Verdict: D (Debated). Similar to the Nevadaplano, this is a debated regional hypothesis.
s4 = "D"

# Statement (5): Laramide-style structures are generally inboard of Sevier-style structures.
# Verdict: C (Consensus). This fundamental spatial relationship is a key observation in Cordilleran geology.
s5 = "C"

# Statement (6): The Sierra Nevada and Idaho Batholiths were formed by the ancestral Cordilleran arc.
# Verdict: C (Consensus). This is a foundational concept of Mesozoic tectonics in North America.
s6 = "C"

# Statement (7): A series of large-magnitude ignimbrite eruptions began in Idaho/Montana and propagated southwest.
# Verdict: C (Consensus). This accurately describes the well-documented track of the Yellowstone hotspot.
s7 = "C"

# Statement (8): The Colorado Plateau reached modern elevations by the late Cretaceous.
# Verdict: D (Debated). The uplift timing and mechanism of the Colorado Plateau are classic, unresolved geologic problems.
s8 = "D"

# Statement (9): The Cascades arc reached modern elevation by the late Oligocene.
# Verdict: D (Debated). This statement is factually inconsistent; the modern high volcanic edifices are much younger.
s9 = "D"

# Statement (10): Laramide-style deformation was caused by subduction of the Shatsky conjugate.
# Verdict: D (Debated). This is a specific, leading hypothesis for the cause of the Laramide orogeny, but it is not a universally accepted fact.
s10 = "D"

# Combine the results into a single string for the final answer.
final_answer = s1 + s2 + s3 + s4 + s5 + s6 + s7 + s8 + s9 + s10

print(final_answer)
