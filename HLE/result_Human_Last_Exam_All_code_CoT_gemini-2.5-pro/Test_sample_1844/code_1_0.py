#
# Description:
# This script calculates the total number of edges on a piece of paper
# after it has been folded four times, its corners cut, and then unfolded.
# The logic involves categorizing the cuts, calculating the number of new edges 
# created by them, and adding the number of segments the original paper's
# edges are broken into.
#

# Part 1: Calculating the number of new edges created directly by the cuts ("cut-edges").
# The four cuts on the folded square correspond to different types of features on the unfolded paper.

# The cut at the 'central' corner (C_c) is mirrored across multiple fold lines,
# creating 4 separate 4-sided holes in the interior of the paper.
holes_from_central_cut = 4
edges_per_hole = 4
edges_from_central_cut = holes_from_central_cut * edges_per_hole

# The cut at the 'original outer' corner (C_oo) snips the four corners of the unfolded paper.
# Each snip creates one new edge on the outer boundary.
snips_from_outer_cut = 4
edges_per_snip = 1
edges_from_outer_cut = snips_from_outer_cut * edges_per_snip

# The two cuts at the 'original edge' corners (C_eo and C_oe) create V-shaped notches
# on the outer boundary of the unfolded paper. Each of these two cuts creates 4 notches,
# and each notch consists of 2 edges.
num_edge_cuts = 2
notches_per_edge_cut = 4
edges_per_notch = 2
edges_from_edge_cuts = num_edge_cuts * notches_per_edge_cut * edges_per_notch

# Sum of all new edges created by the cuts.
total_cut_edges = edges_from_central_cut + edges_from_outer_cut + edges_from_edge_cuts

# Part 2: Calculating the number of edges from the original paper's boundary.
# Each of the 4 original edges is broken into smaller segments by the corner snips and notches.

# Each original edge is interrupted by 2 corner snips (one at each end) and 
# 2 notches (from the two 'edge' cuts). This makes 4 breaks in total.
breaks_per_original_edge = 2 + 2

# The number of segments an edge is broken into is the number of breaks + 1.
segments_per_original_edge = breaks_per_original_edge + 1

# Total segments from the original boundary of the paper.
num_original_edges = 4
total_original_edge_segments = num_original_edges * segments_per_original_edge

# Part 3: Final Calculation.
# The total number of edges is the sum of the new cut-edges and the segments from the original boundary.
total_edges = total_cut_edges + total_original_edge_segments

# Output the detailed calculation.
print("The calculation for the total number of edges is broken down as follows:")
print("\nPart 1: Calculating the number of new edges created directly by the cuts.")
print(f"The cut at the 'central' corner creates {holes_from_central_cut} holes, each with {edges_per_hole} edges.")
print(f"Edges from central cut = {holes_from_central_cut} * {edges_per_hole} = {edges_from_central_cut}")
print(f"The cut at the 'original outer' corner creates {snips_from_outer_cut} snips, each being {edges_per_snip} edge.")
print(f"Edges from outer corner cut = {snips_from_outer_cut} * {edges_per_snip} = {edges_from_outer_cut}")
print(f"The {num_edge_cuts} cuts at the 'original edge' corners each create {notches_per_edge_cut} notches with {edges_per_notch} edges each.")
print(f"Edges from the two edge cuts = {num_edge_cuts} * {notches_per_edge_cut} * {edges_per_notch} = {edges_from_edge_cuts}")
print(f"Total new edges from cuts = {edges_from_central_cut} + {edges_from_outer_cut} + {edges_from_edge_cuts} = {total_cut_edges}")

print("\nPart 2: Calculating the number of edges from the original paper's boundary.")
print("Each of the 4 original edges is interrupted by 2 corner snips and 2 notches, making 4 breaks.")
print(f"Number of segments per original edge = {breaks_per_original_edge} breaks + 1 = {segments_per_original_edge}")
print(f"Total segments from original boundary = {num_original_edges} edges * {segments_per_original_edge} segments/edge = {total_original_edge_segments}")

print("\nPart 3: Final Calculation.")
print("Total edges = (Total new edges from cuts) + (Total segments from original boundary)")
print(f"Total edges = {total_cut_edges} + {total_original_edge_segments} = {total_edges}")