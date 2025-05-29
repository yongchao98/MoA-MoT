def generate_output(input_list):
    output = []
    for i in range(1, len(input_list)):
        if not output or input_list[i] > output[-1]:
            output.append(input_list[i])
    return output

input_list = [20, 32, 57]
output = generate_output(input_list)
print(output)