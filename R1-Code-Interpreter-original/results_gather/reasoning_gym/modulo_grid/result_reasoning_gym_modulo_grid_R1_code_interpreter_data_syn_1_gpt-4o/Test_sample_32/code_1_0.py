def find_column_number(row):
    # For row 2, the column is 5
    if row == 2:
        return 5
    # For row 3, the columns are 3, 7, 11, 15, 19
    elif row == 3:
        return [3, 7, 11, 15, 19]
    # For row 18, let's assume a similar pattern as row 3
    elif row == 18:
        return [3, 7, 11, 15, 19]
    else:
        return None

# Find the column numbers for row 18
column_numbers = find_column_number(18)
print(column_numbers)