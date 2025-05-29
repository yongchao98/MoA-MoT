# Number of classes per weekday and on Saturday
classes_per_weekday = 8387270
classes_on_saturday = 8

# Number of students per class
students_per_class = 18387270

# Charge per student
charge_per_student = 18387270.00

# Total classes in a week
total_classes_week = (classes_per_weekday * 5) + classes_on_saturday

# Total students in a week
total_students_week = total_classes_week * students_per_class

# Total revenue
total_revenue = total_students_week * charge_per_student

print(total_revenue)