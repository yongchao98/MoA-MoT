# Number of separate patches of surface code
num_patches = 2

# The problem states each patch has two holes.
# These are interpreted as internal boundaries.
num_internal_holes_per_patch = 2

# Each planar patch also has one outer boundary.
# The total number of boundaries 'd' is the sum of internal and outer boundaries.
num_boundaries_per_patch = num_internal_holes_per_patch + 1

# The number of logical qubits 'k' for a single patch is d - 1.
logical_qubits_per_patch = num_boundaries_per_patch - 1

# The total number of logical qubits is the sum from each independent patch.
total_logical_qubits = logical_qubits_per_patch * num_patches

print(f"A single patch with {num_internal_holes_per_patch} internal holes has {num_boundaries_per_patch} total boundaries, encoding {logical_qubits_per_patch} logical qubits.")
print(f"For {num_patches} such patches, the total number of logical qubits is:")
print(f"{logical_qubits_per_patch} + {logical_qubits_per_patch} = {total_logical_qubits}")