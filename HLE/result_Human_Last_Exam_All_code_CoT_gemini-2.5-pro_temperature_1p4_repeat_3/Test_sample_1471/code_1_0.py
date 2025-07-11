# The problem involves solving an equation with cuneiform numerals.
# Let's first translate these numerals into modern numbers.
# പ്പെ = 60
# 𒌋 = 10
# 𒁹 = 1

# The first term is പ്പെപ്പെ / 𒌋𒌋𒌋
# Numerator: പ്പെപ്പെ = 60 + 60 = 120
numerator = 120
# Denominator: 𒌋𒌋𒌋 = 10 + 10 + 10 = 30
denominator = 30

# The second term is 𒌋𒌋𒁹𒁹𒁹
# 𒌋𒌋 = 10 + 10 = 20
# 𒁹𒁹𒁹 = 1 + 1 + 1 = 3
# So, the second number is 20 + 3 = 23
second_term = 23

# Now, we calculate the result of the full equation: 120 / 30 - 23
result = (numerator / denominator) - second_term

# Finally, we print the equation in modern numbers along with its solution.
print(f"{numerator} / {denominator} - {second_term} = {int(result)}")