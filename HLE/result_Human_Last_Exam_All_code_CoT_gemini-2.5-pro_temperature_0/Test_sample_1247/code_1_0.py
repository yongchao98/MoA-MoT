# Plan:
# 1. Calculate the number of 1324-avoiding permutations with 3 inversions
#    whose support is a single contiguous block.
# 2. Calculate the number of 1324-avoiding permutations with 3 inversions
#    whose support is split into two disjoint blocks at the ends.
# 3. Sum the results.

# Case 1: Contiguous support
# Based on core permutation '321' (on 3 elements): 2 permutations (at beginning or end)
count_321_like = 2
# Based on core permutation '1432' (on 4 elements): 1 permutation (at the end)
count_1432_like = 1
# Based on core permutation '4123' (on 4 elements): 2 permutations (at beginning or end)
count_4123_like = 2

total_contiguous = count_321_like + count_1432_like + count_4123_like

# Case 2: Disjoint support at the ends (1-inversion block and 2-inversion block)
# Swap at beginning, 2-inv block at end: 2 permutations
count_swap_first = 2
# 2-inv block at beginning, swap at end: 2 permutations
count_swap_last = 2

total_disjoint = count_swap_first + count_swap_last

# Total count
total_av_333_3 = total_contiguous + total_disjoint

print("Calculation for av_{333}^3(1324):")
print("Number of permutations with contiguous support:")
print(f"  - Based on core '321': {count_321_like}")
print(f"  - Based on core '1432': {count_1432_like}")
print(f"  - Based on core '4123': {count_4123_like}")
print(f"Subtotal for contiguous support: {count_321_like} + {count_1432_like} + {count_4123_like} = {total_contiguous}")
print("")
print("Number of permutations with disjoint support at the ends:")
print(f"  - Swap at beginning, 2-inversion block at end: {count_swap_first}")
print(f"  - 2-inversion block at beginning, swap at end: {count_swap_last}")
print(f"Subtotal for disjoint support: {count_swap_first} + {count_swap_last} = {total_disjoint}")
print("")
print("Total number of such permutations:")
print(f"Final equation: {total_contiguous} + {total_disjoint} = {total_av_333_3}")