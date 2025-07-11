m = True
a = set((lambda f: [f.append(f[-1] + f[~m]) or f[-1] for _ in range(len("1 2 3 4 5"[:5]))] and f)([not m, m]))
b = set((lambda f: [f.append(f[~(m & ~m) & m] + f[-1] + f[~m]) or f[-1] for _ in range(len("1 2 3 4 5"[:5]))] and f)([not m]*(m | m)+[m]))
result_set = a - b
print(f"{' + '.join(map(str, sorted(list(result_set))))} = {sum(result_set)}")
