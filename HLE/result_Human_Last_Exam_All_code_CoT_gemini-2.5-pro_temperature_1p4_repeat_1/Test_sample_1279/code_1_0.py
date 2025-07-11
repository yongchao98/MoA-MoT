import math

# This script calculates the symmetry factors for the second-order vacuum bubble
# diagrams in phi^4 theory and sums them.

# --- Diagram 1: Figure-8 Diagram ---

# Description: Two vertices, each with a self-loop, connected by two propagators.
# Symmetries:
# 1. Swapping the two identical vertices.
vertex_swap_sym_fig8 = 2
# 2. Swapping the two identical propagators connecting the vertices.
line_swap_sym_fig8 = 2
# 3. Swapping the ends of the propagator for each of the two loops.
loop_sym_fig8 = 2 * 2

# The total symmetry factor for the Figure-8 diagram is the product of these symmetries.
s1 = vertex_swap_sym_fig8 * line_swap_sym_fig8 * loop_sym_fig8
print("Calculation for the Figure-8 diagram:")
print(f"S_figure_8 = (vertex swaps) * (line swaps) * (loop symmetries)")
print(f"S_figure_8 = {vertex_swap_sym_fig8} * {line_swap_sym_fig8} * {loop_sym_fig8 / 2} * {loop_sym_fig8 / 2} = {s1}")
print("-" * 20)


# --- Diagram 2: Sunset Diagram ---

# Description: Two vertices connected by four propagators.
# Symmetries:
# 1. Swapping the two identical vertices.
vertex_swap_sym_sunset = 2
# 2. Permuting the four identical propagators connecting the vertices.
line_perm_sym_sunset = math.factorial(4)

# The total symmetry factor for the Sunset diagram is the product of these symmetries.
s2 = vertex_swap_sym_sunset * line_perm_sym_sunset
print("Calculation for the Sunset diagram:")
print(f"S_sunset = (vertex swaps) * (permutations of lines)")
print(f"S_sunset = {vertex_swap_sym_sunset} * {line_perm_sym_sunset} = {s2}")
print("-" * 20)


# --- Sum of all symmetry factors ---

total_s = s1 + s2
print("Total sum of symmetry factors:")
print(f"Total S = S_figure_8 + S_sunset")
print(f"Total S = {s1} + {s2} = {total_s}")
