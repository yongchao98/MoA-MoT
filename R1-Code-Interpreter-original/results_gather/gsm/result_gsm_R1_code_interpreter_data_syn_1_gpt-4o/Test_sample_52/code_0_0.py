# Jon's times
jon_swim = 40
jon_bike = 80
jon_run = 50
jon_total = jon_swim + jon_bike + jon_run

# James's times
james_swim = jon_swim * 0.9
james_bike = jon_bike + 5
james_total = jon_total + 10

# Calculate James's run time
james_run = james_total - (james_swim + james_bike)
print(james_run)