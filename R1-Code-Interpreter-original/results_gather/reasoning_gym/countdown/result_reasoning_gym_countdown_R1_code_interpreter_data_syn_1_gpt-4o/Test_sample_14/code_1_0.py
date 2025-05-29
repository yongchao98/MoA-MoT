# Given numbers
numbers = [72, 24, 48, 15]

# Target number
target = 624

# Check multiplication possibilities
for i in range(len(numbers)):
    for j in range(i+1, len(numbers)):
        if numbers[i] * numbers[j] == target:
            print(f"{numbers[i]} * {numbers[j]}")
        if numbers[i] * numbers[j] < target:
            remaining = target - numbers[i] * numbers[j]
            if remaining in numbers:
                print(f"{numbers[i]} * {numbers[j]} + {remaining}")

# Check division possibilities
for number in numbers:
    if target % number == 0:
        quotient = target // number
        if quotient in numbers:
            print(f"{number} * {quotient}")