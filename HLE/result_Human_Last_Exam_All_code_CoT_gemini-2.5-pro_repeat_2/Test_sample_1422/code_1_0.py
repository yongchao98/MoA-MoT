#
# This script calculates the minimum number of vertices required for a two-loop
# Feynman diagram in an interacting scalar field theory.
#

# Step 1: State the fundamental formulas for Feynman diagram topology.
print("To solve this, we use two fundamental topological relations for a connected Feynman diagram:")
print("1. L = I - V + 1")
print("2. n * V = 2 * I + E")
print("\nWhere:")
print("  L = number of loops")
print("  I = number of internal lines (propagators)")
print("  V = number of vertices")
print("  n = number of lines meeting at each vertex (e.g., n=4 for a phi^4 interaction)")
print("  E = number of external lines (incoming/outgoing particles)")
print("-" * 50)

# Step 2: Apply the given conditions.
L = 2
print(f"We are considering a two-loop diagram, so we set L = {L}.")
print("-" * 50)


# Step 3: Combine the formulas to solve for V.
print("First, we rearrange formula (1) to express the number of internal lines (I) in terms of V:")
print("I = L + V - 1")
print(f"Substituting L = {L}, we get: I = {L} + V - 1, which simplifies to I = V + 1")
print("\nNext, we substitute this expression for I into formula (2):")
print("n * V = 2 * (V + 1) + E")
print("n * V = 2*V + 2 + E")
print("\nFinally, we rearrange the equation to solve for the number of vertices (V):")
print("V*n - 2*V = E + 2")
print("V * (n - 2) = E + 2")
print("V = (E + 2) / (n - 2)")
print("-" * 50)

# Step 4: Analyze the formula for V to find its minimum value.
print("To find the minimum possible value for V, we must choose appropriate values for E and n:")
print(" - To minimize V, we must minimize the numerator (E + 2). The smallest non-negative integer for E is 0.")
print(" - To minimize V, we must maximize the denominator (n - 2). For an interacting theory, n must be at least 3.")
print("\nLet's test the simplest interaction theories.")
print("-" * 50)

# Step 5: Test specific cases for n.

# Case 1: phi^3 theory
n_phi3 = 3
E_min = 0
V_num = E_min + 2
V_den_phi3 = n_phi3 - 2
V_phi3 = V_num / V_den_phi3
print(f"Case 1: A phi^{n_phi3} theory (n = {n_phi3}).")
print(f"Using the minimum E = {E_min}, the number of vertices is:")
print(f"V = ({E_min} + 2) / ({n_phi3} - 2) = {V_num} / {V_den_phi3} = {int(V_phi3)}")


# Case 2: phi^4 theory
n_phi4 = 4
V_den_phi4 = n_phi4 - 2
V_phi4 = V_num / V_den_phi4
print(f"\nCase 2: A phi^{n_phi4} theory (n = {n_phi4}).")
print(f"Using the minimum E = {E_min}, the number of vertices is:")
print(f"V = ({E_min} + 2) / ({n_phi4} - 2) = {V_num} / {V_den_phi4} = {int(V_phi4)}")
print("-" * 50)

# Step 6: Conclude with the final answer.
min_V = min(V_phi3, V_phi4)
print("By comparing the cases, we can see that a phi^4 interaction allows for a two-loop diagram with fewer vertices than a phi^3 interaction.")
print(f"\nThe minimum number of vertices is therefore {int(min_V)}.")
print("This corresponds to the 'figure-8' vacuum diagram in a phi^4 theory.")
