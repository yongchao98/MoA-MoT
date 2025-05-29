# Jennifer's circle
jennifer_male_colleagues = 1
jennifer_female_colleagues = 6

# William's circle
william_male_colleagues = 3
william_female_colleagues = 3

# Matilda is part of both circles
# Excluding Matilda from Jennifer's circle
jennifer_female_colleagues_excluding_matilda = jennifer_female_colleagues - 1

# Excluding Matilda from William's circle
william_female_colleagues_excluding_matilda = william_female_colleagues - 1

# Total unique female colleagues for Matilda
matilda_female_colleagues = jennifer_female_colleagues_excluding_matilda + william_female_colleagues_excluding_matilda

print(matilda_female_colleagues)