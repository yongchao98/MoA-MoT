# The number of loops in the Feynman diagram.
L = 2

# To find the minimum number of vertices, we must consider the diagram
# with the minimum possible number of external lines (E).
# A vacuum diagram has E=0, which is the minimum possible value.
E_min = 0

# For a scalar phi^4 theory, the number of vertices (V) is related to
# the number of loops (L) and external lines (E) by the formula:
# V = L - 1 + E / 2

# We plug in the values L=2 and E=0 to find the minimum number of vertices.
V_min = L - 1 + E_min / 2

print("The formula for the number of vertices (V) in a phi^4 theory is: V = L - 1 + E / 2")
print("Where L is the number of loops and E is the number of external lines.")
print("\nWe are looking for a two-loop diagram, so L = 2.")
print("To find the *minimum* number of vertices, we use the minimum number of external lines, E = 0.")
print("\nPlugging these values into the formula:")
# The final equation with each number explicitly shown
print(f"V_min = {L} - 1 + {int(E_min)} / 2")
print(f"V_min = {int(V_min)}")

print("\nTherefore, the minimum number of vertices is 1.")
print("This corresponds to a 'figure-eight' vacuum diagram where two loops are attached to a single vertex.")
