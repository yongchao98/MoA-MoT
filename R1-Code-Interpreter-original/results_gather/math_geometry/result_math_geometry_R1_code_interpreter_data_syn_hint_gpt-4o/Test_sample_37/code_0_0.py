count = 0
for a in range(6, 12):
    b = 23 - 2 * a
    if b > 0 and 2 * a > b:
        count += 1
print(count)