# Jessica's circle
jessica_female_colleagues = 1  # excluding Matilda

# Michael's circle
michael_female_colleagues = 6  # including Matilda

# Matilda's female colleagues
# Matilda is part of both circles, so we need to consider the unique female colleagues
# in Michael's circle excluding Matilda herself
matilda_female_colleagues = jessica_female_colleagues + (michael_female_colleagues - 1)

print(matilda_female_colleagues)