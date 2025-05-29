test_input = [5, 0, 5, 5, 5, 0, 5, 5, 5, 5, 0, 5, 5, 0, 0, 0, 0]
test_output = test_input[:]
test_output[1] = 0
test_output[4] = 0
print("<<<" + " ".join(map(str, test_output)) + ">>>")