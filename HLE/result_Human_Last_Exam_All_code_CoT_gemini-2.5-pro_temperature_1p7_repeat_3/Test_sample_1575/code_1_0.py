# Main parameters
N = 100  # Number of elements
S = 5    # Non-adjacent swap distance

# Number of elements per group
group_size = N // S

print("This problem can be solved by analyzing the movement of elements between groups defined by their position modulo 5.")
print(f"There are {N} elements and non-adjacent swaps can happen between elements at positions i and i+{S}.")
print("\nStep 1: Group positions and elements.")
print("The 'free' non-adjacent swaps mean we can rearrange elements within a set of positions {k, k+S, k+2S, ...} at no cost.")
print(f"This partitions the {N} positions into {S} groups, P_0, P_1, ..., P_{S-1}, based on position modulo {S}.")
print(f"Each position group P_k has {group_size} positions.\n")
print("Similarly, we can group the elements by their value modulo 5: V_0, V_1, ..., V_4.")
print("V_m contains elements {x | 1 <= x <= 100, x mod 5 = m}. Note V_0 are multiples of 5.")
print(f"Each value group V_m has {group_size} elements.\n")

print("Step 2: Determine initial and final configurations of the position groups.")
print("Initial state: The element at position 'j' is 'j+1'.")
print("  - P_0 (pos 0,5,..) holds elements {1,6,..} -> V_1")
print("  - P_1 (pos 1,6,..) holds elements {2,7,..} -> V_2")
print("  - P_2 (pos 2,7,..) holds elements {3,8,..} -> V_3")
print("  - P_3 (pos 3,8,..) holds elements {4,9,..} -> V_4")
print("  - P_4 (pos 4,9,..) holds elements {5,10,..} -> V_0\n")

print("Final state: The sequence is reversed. The element at position 'j' is '100-j'.")
print("  - P_0 must hold elements {100-j | j = 0 mod 5}. These values are {100,95,..} -> V_0")
print("  - P_1 must hold elements {100-j | j = 1 mod 5}. These values are {99,94,..} -> V_4")
print("  - P_2 must hold elements {100-j | j = 2 mod 5}. These values are {98,93,..} -> V_3")
print("  - P_3 must hold elements {100-j | j = 3 mod 5}. These values are {97,92,..} -> V_2")
print("  - P_4 must hold elements {100-j | j = 4 mod 5}. These values are {96,91,..} -> V_1\n")

print("Step 3: Calculate moves required by counting swaps across boundaries.")
print("The 'adjacent swaps' move elements between adjacent position groups (P_k <-> P_{k+1}).")
print("The position groups form a line: P_0 - P_1 - P_2 - P_3 - P_4.")
print("The total moves is the sum of swaps required at each of the 4 boundaries.\n")

total_swaps = 0
# Boundary P_0 <-> P_1
left_initial_0 = {'V_1'}
left_final_0 = {'V_0'}
moved_out_0 = left_initial_0 - left_final_0
moved_in_0 = left_final_0 - left_initial_0
swaps_01 = len(moved_out_0) * group_size
print(f"Boundary P_0 <-> P_1:")
print(f"  - Initially, P_0 holds {left_initial_0}. Finally, it must hold {left_final_0}.")
print(f"  - This means {moved_out_0} must cross to the right, and {moved_in_0} must cross to the left.")
print(f"  - Number of swaps needed = {swaps_01}\n")
total_swaps += swaps_01

# Boundary P_1 <-> P_2
left_initial_1 = {'V_1', 'V_2'}
left_final_1 = {'V_0', 'V_4'}
moved_out_1 = left_initial_1 - left_final_1
moved_in_1 = left_final_1 - left_initial_1
swaps_12 = len(moved_out_1) * group_size
print(f"Boundary P_1 <-> P_2:")
print(f"  - Initially, {{P_0,P_1}} hold {left_initial_1}. Finally, they must hold {left_final_1}.")
print(f"  - This means {moved_out_1} must cross to the right, and {moved_in_1} must cross to the left.")
print(f"  - Number of swaps needed = {swaps_12}\n")
total_swaps += swaps_12

# Boundary P_2 <-> P_3
left_initial_2 = {'V_1', 'V_2', 'V_3'}
left_final_2 = {'V_0', 'V_4', 'V_3'}
moved_out_2 = left_initial_2 - left_final_2
moved_in_2 = left_final_2 - left_initial_2
swaps_23 = len(moved_out_2) * group_size
print(f"Boundary P_2 <-> P_3:")
print(f"  - Initially, {{P_0,P_1,P_2}} hold {left_initial_2}. Finally, they must hold {left_final_2}.")
print(f"  - This means {moved_out_2} must cross to the right, and {moved_in_2} must cross to the left.")
print(f"  - Number of swaps needed = {swaps_23}\n")
total_swaps += swaps_23

# Boundary P_3 <-> P_4
left_initial_3 = {'V_1', 'V_2', 'V_3', 'V_4'}
left_final_3 = {'V_0', 'V_4', 'V_3', 'V_2'}
moved_out_3 = left_initial_3 - left_final_3
moved_in_3 = left_final_3 - left_initial_3
swaps_34 = len(moved_out_3) * group_size
print(f"Boundary P_3 <-> P_4:")
print(f"  - Initially, {{P_0,P_1,P_2,P_3}} hold {left_initial_3}. Finally, they must hold {left_final_3}.")
print(f"  - This means {moved_out_3} must cross to the right, and {moved_in_3} must cross to the left.")
print(f"  - Number of swaps needed = {swaps_34}\n")
total_swaps += swaps_34


print("Step 4: Sum the swaps for the final answer.")
print(f"Total minimum moves = {swaps_01} + {swaps_12} + {swaps_23} + {swaps_34} = {total_swaps}")