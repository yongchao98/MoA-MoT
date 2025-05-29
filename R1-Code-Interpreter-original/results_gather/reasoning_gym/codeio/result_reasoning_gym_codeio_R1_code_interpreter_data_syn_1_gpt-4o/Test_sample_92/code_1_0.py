def total_plates(pattern):
    total = 1
    for p in pattern:
        if p == "N":
            total *= 10
        elif p == "L":
            total *= 26
    return total

def main_solution(population):
    patterns = {0: [(1, "")]}
    pattern_length = 1
    while True:
        old = patterns[pattern_length - 1]
        new = []
        for (x, p) in old:
            new.append((total_plates(p + "N"), p + "N"))
            new.append((total_plates(p + "L"), p + "L"))
        new = list(dict(new).items())
        patterns[pattern_length] = new
        pattern_length += 1
        if any(t >= population for t, p in new):
            break
    
    tps = []
    for tp in patterns.values():
        tps += [(t, p) for (t, p) in tp]
    patterns = dict(tps)
    sortedkeys = sorted(list(patterns.keys()))
    
    for i, t in enumerate(sortedkeys):
        if t >= population:
            pattern = patterns[t]
            total_plates = t
            excess_plates = t - population
            previous_pattern = patterns[sortedkeys[i - 1]] if i > 0 else ""
            next_pattern = patterns[sortedkeys[i + 1]] if i + 1 < len(sortedkeys) else ""
            return {
                "pattern": pattern,
                "total_plates": total_plates,
                "excess_plates": excess_plates,
                "previous_pattern": previous_pattern,
                "next_pattern": next_pattern
            }

# Example usage
result = main_solution(574888)
print(result)