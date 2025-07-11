# Number of '0'-blocks in the worst-case initial sequence '0101...'.
# To transform this to '1111...', we need to remove all '0'-blocks.
ops_to_transform_initial_to_monochrome = 50

# Number of '0'-blocks in the worst-case target sequence '1010...'.
# To transform '1111...' to this target, we start with the first '1'-block
# and need to insert all the '0'-blocks.
ops_to_transform_monochrome_to_target = 50

# The total minimum number of operations in the worst-case scenario.
n = ops_to_transform_initial_to_monochrome + ops_to_transform_monochrome_to_target

# Output the final equation
print(f"{ops_to_transform_initial_to_monochrome} + {ops_to_transform_monochrome_to_target} = {n}")