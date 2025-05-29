count = 0
for b in range(1, 10):  # b must be between 1 and 9
    for a in range(0, b):  # a must be less than b
        for c in range(0, b):  # c must be less than b
            if a != 0:  # a must be non-zero for a 3-digit number
                count += 1

print(count)