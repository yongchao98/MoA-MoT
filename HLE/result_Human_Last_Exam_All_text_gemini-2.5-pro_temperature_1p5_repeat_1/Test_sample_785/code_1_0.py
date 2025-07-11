count = 0
for n6 in range(1000 // 6 + 1):
    rem_after_6 = 1000 - 6 * n6
    for n5 in range(rem_after_6 // 5 + 1):
        rem_after_5 = rem_after_6 - 5 * n5
        for n4 in range(rem_after_5 // 4 + 1):
            n1 = rem_after_5 - 4 * n4
            count += (n1 + 1) * (n4 + 1) * (n5 + 1)