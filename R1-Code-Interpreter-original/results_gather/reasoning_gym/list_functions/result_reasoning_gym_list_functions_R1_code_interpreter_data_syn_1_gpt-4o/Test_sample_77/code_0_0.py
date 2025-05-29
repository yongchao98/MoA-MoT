def determine_output(input_list):
    increasing = 0
    decreasing = 0
    
    for i in range(1, len(input_list)):
        if input_list[i] > input_list[i - 1]:
            increasing += 1
        elif input_list[i] < input_list[i - 1]:
            decreasing += 1
    
    if increasing > decreasing:
        return [1]
    else:
        return [0]

input_data = [44, 80, 94]
output = determine_output(input_data)
print(output)