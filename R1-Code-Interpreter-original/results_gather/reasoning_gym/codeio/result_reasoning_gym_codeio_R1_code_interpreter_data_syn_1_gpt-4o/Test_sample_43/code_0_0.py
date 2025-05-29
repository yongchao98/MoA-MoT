def generate_demlo_numbers(n):
    result = []
    for i in range(1, n + 1):
        demlo_number = ((10**i) // 9)**2
        result.append(demlo_number)
    return result

# Given input
n = 7
# Generate Demlo numbers
output = generate_demlo_numbers(n)
print(output)