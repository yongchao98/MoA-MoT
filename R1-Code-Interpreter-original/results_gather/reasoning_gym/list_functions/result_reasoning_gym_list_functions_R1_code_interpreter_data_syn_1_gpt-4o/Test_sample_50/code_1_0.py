def find_closest_to_median(input_list):
    sorted_list = sorted(input_list)
    n = len(sorted_list)
    median = (sorted_list[n // 2] if n % 2 == 1 else (sorted_list[n // 2 - 1] + sorted_list[n // 2]) / 2)
    closest_to_median = min(input_list, key=lambda x: abs(x - median))
    return closest_to_median

# Apply the rule to the new input
output_element = find_closest_to_median([69, 48, 2])
print(output_element)