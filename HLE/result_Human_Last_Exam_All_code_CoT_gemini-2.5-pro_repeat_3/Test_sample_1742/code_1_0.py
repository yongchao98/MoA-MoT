import numpy as np

# Dimension vectors for the indecomposable modules of A=C(1->2->3)
P1 = np.array([1, 1, 1])
S2 = np.array([0, 1, 0])
I2 = np.array([1, 1, 0])

# The tau-tilting module in question is T = P1 + S2 + I2
T = P1 + S2 + I2

print("The unique tau-tilting module that is not a slice is T = P1 + S2 + I2.")
print("Its indecomposable summands have dimension vectors:")
print(f"P1: ({P1[0]},{P1[1]},{P1[2]})")
print(f"S2: ({S2[0]},{S2[1]},{S2[2]})")
print(f"I2: ({I2[0]},{I2[1]},{I2[2]})")

# The final answer is the module itself, represented by its summands
print("\nThe module is P1 (+) S2 (+) I2.")
print(f"P1 is represented by ({P1[0]}, {P1[1]}, {P1[2]})")
print(f"S2 is represented by ({S2[0]}, {S2[1]}, {S2[2]})")
print(f"I2 is represented by ({I2[0]}, {I2[1]}, {I2[2]})")
