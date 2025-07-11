# Number of loops in the diagram.
L = 2

# To find the minimum number of vertices (V), we must use the minimum
# possible number of external lines (E). The absolute minimum is E=0, which
# corresponds to a vacuum diagram.
E = 0

# The number of vertices (V) in a phi^4 theory is related to the number
# of loops (L) and external lines (E) by the formula: V = L - 1 + E / 2.
# We will now calculate the minimum number of vertices.
min_V = L - 1 + E / 2

print("For a Feynman diagram in an interacting scalar field theory (phi^4):")
print(f"The number of vertices (V) is given by the formula: V = L - 1 + E / 2")
print(f"Where L is the number of loops and E is the number of external lines.\n")
print("To find the minimum number of vertices for L=2 loops, we must use the minimum number of external lines, which is E=0.")
print("\nCalculating the minimum number of vertices:")
print(f"V_min = {L} - 1 + {E} / 2")
print(f"V_min = {int(min_V)}")
print(f"\nThe minimum number of vertices is {int(min_V)}. This corresponds to the 'figure-eight' vacuum diagram.")
