# Total number of students
total_students = 2780965

# Students through exit A
students_exit_A = 0.30 * total_students

# Remaining students after exit A
remaining_students_after_A = total_students - students_exit_A

# Students through exit B
students_exit_B = (3/5) * remaining_students_after_A

# Students through exit C
students_exit_C = remaining_students_after_A - students_exit_B

# Output the number of students who went out through exit C
print(int(students_exit_C))