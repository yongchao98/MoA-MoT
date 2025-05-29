# Number of students needing each size
extra_small = 5351035
small = 2 * extra_small
medium = small - 4
large = medium // 2
extra_large = large + 6

# Total number of shirts
total_shirts = extra_small + small + medium + large + extra_large

print(total_shirts)