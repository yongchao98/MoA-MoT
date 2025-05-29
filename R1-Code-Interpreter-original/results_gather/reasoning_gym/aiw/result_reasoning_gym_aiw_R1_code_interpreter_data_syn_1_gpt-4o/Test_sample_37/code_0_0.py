# Patricia's circle
patricia_female_colleagues = 5

# Robert's circle
robert_female_colleagues = 3

# Matilda is counted in both circles
# Excluding Matilda from each circle
matilda_female_colleagues = (patricia_female_colleagues - 1) + (robert_female_colleagues - 1)

print(matilda_female_colleagues)