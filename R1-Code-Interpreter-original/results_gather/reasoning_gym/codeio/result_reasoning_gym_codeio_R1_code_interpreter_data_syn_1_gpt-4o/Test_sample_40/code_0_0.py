def sort_pebbles(arr):
    red, white, blue = 0, 0, len(arr)
    
    while white < blue:
        if arr[white] == 'r':
            arr[red], arr[white] = arr[white], arr[red]
            red += 1
            white += 1
        elif arr[white] == 'w':
            white += 1
        else:  # arr[white] == 'b'
            blue -= 1
            arr[white], arr[blue] = arr[blue], arr[white]
    
    return arr

# Given input
arr = ['b', 'r', 'w', 'b', 'b', 'b', 'w', 'w', 'b', 'b', 'b', 'b', 'r', 'b', 'b', 'w', 'r', 'w', 'r', 'b']
sorted_arr = sort_pebbles(arr)
print(sorted_arr)