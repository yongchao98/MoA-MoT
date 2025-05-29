count = 0
for a in range(1, 10):  # tens digit from 1 to 9
    for b in range(0, 10):  # ones digit from 0 to 9
        if a > b:
            count += 1
print(count)