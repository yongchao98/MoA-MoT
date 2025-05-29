def find_elements(input_list):
    # Sort the list to find the largest and smallest elements
    sorted_list = sorted(input_list)
    # Select the largest and smallest elements
    largest_elements = sorted_list[-3:]  # Three largest elements
    smallest_elements = sorted_list[:2]  # Two smallest elements
    # Combine them in descending order
    result = sorted(largest_elements + smallest_elements, reverse=True)
    return result

# Test the function with the new input
input_list = [38, 49, 41, 28, 27, 58, 69]
output = find_elements(input_list)
print(output)