def find_50th_segmented_number():
    segmented_numbers = [1]
    num_to_check = 2
    while len(segmented_numbers) < 50:
        is_sum_of_consecutive = False
        for i in range(len(segmented_numbers)):
            for j in range(i + 2, len(segmented_numbers) + 1):
                current_sum = sum(segmented_numbers[i:j])
                if current_sum == num_to_check:
                    is_sum_of_consecutive = True
                    break
                if current_sum > num_to_check:
                    break
            if is_sum_of_consecutive:
                break
        
        if not is_sum_of_consecutive:
            segmented_numbers.append(num_to_check)
            
        num_to_check += 1
        
    return segmented_numbers[49]

fiftyth_element = find_50th_segmented_number()
print(fiftyth_element)