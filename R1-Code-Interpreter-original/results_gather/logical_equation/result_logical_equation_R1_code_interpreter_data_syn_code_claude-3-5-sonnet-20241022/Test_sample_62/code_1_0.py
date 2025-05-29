def find_solution():
    # Given numbers
    numbers = {1, 2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50}
    result = [0] * 13  # A to M
    
    # From G = 3.6J
    # Since J must be a whole number and G must be in our set
    # Let's find possible J and G pairs
    for j in numbers:
        g = 3.6 * j
        if g in numbers:
            # We found J and G
            result[9] = j  # J
            result[6] = int(g)  # G
            numbers.remove(j)
            numbers.remove(g)
            break
    
    # From A - J = 6
    # Since we know J, we can find A
    a = result[9] + 6
    if a in numbers:
        result[0] = a  # A
        numbers.remove(a)
    
    # From B - A = -9
    # Since we know A, we can find B
    b = result[0] - 9
    if b in numbers:
        result[1] = b  # B
        numbers.remove(b)
    
    # From B + C = 52
    # Since we know B, we can find C
    c = 52 - result[1]
    if c in numbers:
        result[2] = c  # C
        numbers.remove(c)
    
    # From C = 3.0E
    # Since we know C, we can find E
    e = c / 3.0
    if e in numbers:
        result[4] = int(e)  # E
        numbers.remove(e)
    
    # From E + I = 18
    # Since we know E, we can find I
    i = 18 - result[4]
    if i in numbers:
        result[8] = i  # I
        numbers.remove(i)
    
    # From D - E = -13
    # Since we know E, we can find D
    d = result[4] - 13
    if d in numbers:
        result[3] = d  # D
        numbers.remove(d)
    
    # From K - C = 5
    # Since we know C, we can find K
    k = result[2] + 5
    if k in numbers:
        result[10] = k  # K
        numbers.remove(k)
    
    # From L - B = 17
    # Since we know B, we can find L
    l = result[1] + 17
    if l in numbers:
        result[11] = l  # L
        numbers.remove(l)
    
    # From B + H = 12
    # Since we know B, we can find H
    h = 12 - result[1]
    if h in numbers:
        result[7] = h  # H
        numbers.remove(h)
    
    # The remaining number must be F
    result[5] = numbers.pop()  # F
    
    # The last remaining number must be M
    result[12] = numbers.pop()  # M
    
    print(result)

find_solution()