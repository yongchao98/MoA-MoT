# Number of classes per weekday and on Saturday
classes_per_weekday = 8387270
classes_on_saturday = 8

# Total number of classes in a week
total_classes_weekdays = classes_per_weekday * 5
total_classes_week = total_classes_weekdays + classes_on_saturday

# Number of students per class
students_per_class = 18387270

# Total number of students in a week
total_students_week = total_classes_week * students_per_class

# Charge per student
charge_per_student = 18387270.00

# Total earnings in a week
total_earnings_week = total_students_week * charge_per_student

print(total_earnings_week)