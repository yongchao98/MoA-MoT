import math

def find_rows_with_value(value):
    rows = set()
    n = 0
    while True:
        found = False
        for k in range(n + 1):
            if math.comb(n, k) == value:
                rows.add(n)
                found = True
        if not found and n > value:
            break
        n += 1
    return rows

rows_with_43 = find_rows_with_value(43)
print(len(rows_with_43))