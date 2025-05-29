# Number of students needing extra-small shirts
extra_small = 5351035

# Number of students needing small shirts
small = 2 * extra_small

# Number of students needing medium shirts
medium = small - 4

# Number of students needing large shirts
large = medium // 2

# Number of students needing extra-large shirts
extra_large = large + 6

# Total number of shirts ordered
total_shirts = extra_small + small + medium + large + extra_large

print(total_shirts)