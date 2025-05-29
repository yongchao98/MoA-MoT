def calculate_distinct_terms(n):
    collect_powers = set()
    
    for a in range(2, n + 1):
        for b in range(2, n + 1):
            current_pow = a ** b
            collect_powers.add(current_pow)
    
    return {"distinct_terms": len(collect_powers)}

# Given input
n = 29
# Calculate the result
result = calculate_distinct_terms(n)
# Print the result
print(result)