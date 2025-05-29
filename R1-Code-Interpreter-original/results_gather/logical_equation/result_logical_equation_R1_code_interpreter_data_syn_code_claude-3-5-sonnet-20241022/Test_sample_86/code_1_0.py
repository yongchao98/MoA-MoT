def find_solution():
    numbers = set([2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96])
    
    # From G + I = 13, both must be small numbers
    for i in numbers:
        for g in numbers:
            if i + g == 13:
                # From J - I = 4
                j = i + 4
                if j in numbers:
                    # From A = 4.0J
                    a = int(4.0 * j)
                    if a in numbers:
                        # From F + G = 60
                        for f in numbers:
                            if f + g == 60:
                                # From F + J = 57
                                if f + j == 57:
                                    # From F + L = 95
                                    l = 95 - f
                                    if l in numbers:
                                        # From J = 1.4K
                                        k = round(j/1.4)
                                        if k in numbers:
                                            # From H = 2.4E
                                            for e in numbers:
                                                h = int(2.4 * e)
                                                if h in numbers:
                                                    # From C + H = 52
                                                    c = 52 - h
                                                    if c in numbers:
                                                        # From D + G = 34
                                                        d = 34 - g
                                                        if d in numbers:
                                                            remaining = numbers - {a,c,d,e,f,g,h,i,j,k,l}
                                                            if len(remaining) == 2:
                                                                b, m = sorted(list(remaining))
                                                                result = [a,b,c,d,e,f,g,h,i,j,k,l,m]
                                                                print(result)
                                                                return

find_solution()