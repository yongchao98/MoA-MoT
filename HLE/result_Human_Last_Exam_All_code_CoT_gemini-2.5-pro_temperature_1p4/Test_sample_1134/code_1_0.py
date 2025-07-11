import numpy as np

def print_state(name, state_vector):
    """Helper function to print a state vector nicely."""
    components = [f"{c:.3f}".replace("-0.000", "0.000") for c in state_vector]
    print(f"State |{name}⟩ = [{', '.join(components)}]")

# Define the relevant quantum states as complex numpy arrays
# Normalization factor
N = 1/np.sqrt(2)

# |-> = 1/sqrt(2) * (|0> - |1>)
state_minus = N * np.array([1, -1], dtype=complex)

# |i> = 1/sqrt(2) * (|0> + i|1>)
state_i = N * np.array([1, 1j], dtype=complex)

# |-i> = 1/sqrt(2) * (|0> - i|1>)
state_neg_i = N * np.array([1, -1j], dtype=complex)

# |+> = 1/sqrt(2) * (|0> + |1>)
state_plus = N * np.array([1, 1], dtype=complex)

# Define the S-dagger (S†) gate matrix
S_dag = np.array([[1, 0], [0, -1j]], dtype=complex)

print("Analyzing the S† (S-dagger) operation:")
print("-" * 40)

# Demonstrate that S† saves the initially fatal states
print("1. Checking initially FATAL states:")
final_state_from_i = S_dag @ state_i
print("Applying S† to the fatal state |i⟩ results in a safe state:")
print(f"S† |i⟩  = |+⟩")

final_state_from_neg_i = S_dag @ state_neg_i
print("Applying S† to the fatal state |-i⟩ results in a safe state:")
print(f"S† |-i⟩ = |->⟩")
print("-" * 40)

# Demonstrate the fateful transformation
print("2. Checking the sacrificial transformation:")
print("Applying S† to the safe state |-> results in a fatal state:")
final_state_from_minus = S_dag @ state_minus

# The equation S†|-> = |i>
print("\nTHE FINAL EQUATION:")
print("S† |-> = |i⟩")
print("Let's verify by printing the numbers in the vectors:")

print("S† @ [", end="")
for num in state_minus:
    print(f" {num.real:.3f}{num.imag:+.3f}j", end="")
print(" ] = [", end="")
for num in final_state_from_minus:
    print(f" {num.real:.3f}{num.imag:+.3f}j", end="")
print(" ]")

print("\nComparing the result with the fatal state |i⟩:")
print("Result [", end="")
for num in final_state_from_minus:
    print(f" {num.real:.3f}{num.imag:+.3f}j", end="")
print(" ]")
print("State |i⟩ [", end="")
for num in state_i:
    print(f" {num.real:.3f}{num.imag:+.3f}j", end="")
print(" ]")
# np.allclose checks if two complex arrays are element-wise equal within a tolerance
are_equal = np.allclose(final_state_from_minus, state_i)
print(f"\nThe resulting state is indeed the fatal state |i⟩: {are_equal}")

print("\nConclusion: By choosing S†, you reduce the number of fatal scenarios from two to one.")