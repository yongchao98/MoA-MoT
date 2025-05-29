def total_plates(pattern):
    total = 1
    for p in pattern:
        if p == "N":
            total *= 10
        elif p == "L":
            total *= 26
    return total

def main_solution(population):
    # Start with the simplest pattern
    patterns = [""]
    while True:
        new_patterns = []
        for pattern in patterns:
            new_patterns.append(pattern + "N")
            new_patterns.append(pattern + "L")
        
        # Calculate total plates for each new pattern
        pattern_totals = [(total_plates(p), p) for p in new_patterns]
        pattern_totals.sort()
        
        # Find the simplest pattern that meets or exceeds the population
        for total, pattern in pattern_totals:
            if total >= population:
                # Find previous and next patterns
                previous_pattern = ""
                next_pattern = ""
                for t, p in pattern_totals:
                    if t < population:
                        previous_pattern = p
                    if t > total:
                        next_pattern = p
                        break
                
                return {
                    "pattern": pattern,
                    "total_plates": total,
                    "excess_plates": total - population,
                    "previous_pattern": previous_pattern,
                    "next_pattern": next_pattern
                }
        
        # Update patterns for the next iteration
        patterns = new_patterns

# Given input
population = 574888
result = main_solution(population)
print(result)